options(warn=2)

suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pbapply))

pboptions(type="timer")

# Thanks to https://www.r-bloggers.com/much-faster-unnesting-with-data-table/
unnest_dt <- function(tbl, col) {
  tbl <- as.data.table(tbl)
  col <- ensyms(col)
  clnms <- syms(setdiff(colnames(tbl), as.character(col)))
  tbl <- as.data.table(tbl)
  tbl <- eval(
    expr(tbl[, as.character(unlist(!!!col)), by = list(!!!clnms)])
  )
  colnames(tbl) <- c(as.character(clnms), as.character(col))
  tbl
}


load(snakemake@input$tax)
 infiles <- snakemake@input$data
 proteomefile <- snakemake@input$proteomes
 alignmentfiles <- snakemake@input$aln
 threads <- snakemake@threads
 domain_compl <- snakemake@params$domaincompl
 marker_compl <- snakemake@params$markercompl
 
 names(alignmentfiles) <- snakemake@params$names
 names(infiles) <- snakemake@params$names
 setDTthreads(1)
 taxcols = c("superkingdom",
             "phylum",
             "class",
             "order",
             "family",
             "genus",
             "species")


 readProteins <- function(file){
   marker_seqdata <- fread(file, select=c("Entry","Proteomes"),index="Entry")
   marker_seqdata[,or:=Proteomes] 
   marker_seqdata[,Proteomes:=strsplit(Proteomes,",",fixed=T)]
   marker_seqdata <- unnest_dt(marker_seqdata, "Proteomes") 
   marker_seqdata[,Proteomes:=trimws(Proteomes)]
   marker_seqdata[,Proteomes:=sub("(^[^:]+): (.+)$","\\1~\\2",Proteomes)]
   marker_seqdata[,(c("Proteome","Chromosome")):=tstrsplit(Proteomes,split="~",fixed=T)]
   marker_seqdata[,Proteomes:=NULL]
   marker_seqdata
 }
 
 marker_seqdata <- rbindlist(pblapply(names(infiles),function(marker)cbind(marker,readProteins(infiles[marker])),cl=threads))
 
 marker_stats <- marker_seqdata[, .(entrycount=.N, proteomecount=length(unique(Proteome))), by=marker]
 marker_multicounts <- as.data.table(dcast(marker_seqdata, "Proteome~marker", value.var="Entry", fun.aggregate=function(x)length(x)>1, fill=0)[,sapply(.SD,sum),.SD=names(infiles)],keep.rownames=T)
 colnames(marker_multicounts) <- c("marker", "multicounts")
 marker_stats <- merge(marker_multicounts, marker_stats)
 
 proteomes <- dcast(marker_seqdata, "Proteome~marker", value.var="Entry", fun.aggregate=list, fill=list())
 proteomes_meta <- fread(proteomefile)
 setnames(proteomes_meta,"Proteome ID","Proteome")
 proteomes <- merge(proteomes, proteomes_meta,by="Proteome")
 rm(proteomes_meta)
 proteomes[,Organism:=NULL]
 proteomes <- cbind(proteomes, lookup(proteomes$`Organism ID`,taxcols,threads=threads))
 rm(tax_graph,tax_names)
 proteomes[,query:=NULL]
 setcolorder(proteomes, c("Proteome",names(infiles)))
 
 superkingdomcompl <- proteomes[,as.data.table(lapply(.SD, function(col)sum(sapply(col,length)!=0)))/.N*100,.SD=names(infiles),by=superkingdom]
 selmarker <- data.table(marker=names(infiles),included=superkingdomcompl[,sapply(.SD, function(col)all(col>domain_compl)),.SD=names(infiles)])
 markercompl <- dcast(melt(superkingdomcompl,id="superkingdom",variable.name="marker"),marker ~ superkingdom)
 marker_stats <- merge(marker_stats, selmarker)
 marker_stats <- merge(marker_stats, markercompl)
 marker_stats[ , marker:=factor(marker,levels=names(infiles),ordered=T)]
 setorder(marker_stats, marker)
 selmarker <- marker_stats[included==T, as.character(marker)]
 proteomes[,multicopy_perc:=100*rowSums(as.data.frame(lapply(.SD, function(col)sapply(col,length)>1)))/length(.SD),.SD=selmarker]
 proteomes[,completion_perc:=100*rowSums(as.data.frame(lapply(.SD, function(col)sapply(col,function(el)length(el)>0))))/length(.SD),.SD=selmarker]
 proteomes[,included:=completion_perc>marker_compl]
 
 
 entries <- marker_seqdata[,.(list(Entry)), by=marker]
 rm(marker_seqdata)
 overlap <- rbindlist(lapply(combn(nrow(entries), m=2, simplify=F), function(x)data.table(entries[x[[1]],marker],entries[x[[2]],marker],overlapping=list(unlist(intersect(entries[x[[1]],unlist(V1)], entries[x[[2]],unlist(V1)]))))))
 overlap[,overlapn:=sapply(overlapping,length)]
 overlapmat <- matrix(NA, nrow=length(infiles), ncol=length(infiles), dimnames=list(names(infiles), names(infiles)))
 overlapmat[as.matrix(overlap[, 1:2])] <- overlap$overlapn
 overlapmat[lower.tri(overlapmat)] <- t(overlapmat)[lower.tri(overlapmat)]



JC <- function(seq1vec,seq2vec){
   gapcorrection <- seq1vec != "-" | seq2vec != "-"
   identitcalpos <- seq1vec[gapcorrection] == seq2vec[gapcorrection]
   identity <- sum(identitcalpos)/sum(gapcorrection)
   b <- 19/20
   if(identity <= (1-b)){
      return(Inf)
   }
   jc <- -b * log(1-(1-identity)/b)
   jc
}
 
 JCs <- function(seqlistrows,seqlistcols){
    grid <- expand.grid(row=seq_along(seqlistrows),col=seq_along(seqlistcols))
    grid$dist <- apply(grid,1,function(i)JC(seqlistrows[[i[[1]]]],seqlistcols[[i[[2]]]]))
    mat <- matrix(grid$dist,nrow=length(seqlistrows),ncol=length(seqlistcols),byrow=F)
    mat
 }
 
 set.seed(123)
 
# From a list of multicandidates list(c(A,B,C),c(E,D)) -> c(B,E) with B and E being the best candidates
selectCandidate <-  function(entries,alignmentfile,threads){
  alignment.look <- read.alignment(alignmentfile,"fasta")
  alignment.look <- as.data.table(unclass(alignment.look)[c("nam","seq")],key="nam")
  alignment.look[,seq:=unlist(seq)]
  overlapping <- alignment.look[unique(unlist(entries))]
  overlapping[,seq:=list(lapply(seq,s2c))]
  sampled <- alignment.look[sample(nam,min(1000, length(nam))), .(nam,seq=lapply(seq,s2c))]
  calc <- function(seqA){
    distvals <- sapply(sampled$seq,function(seqB)JC(seqA,seqB))
    data.table(mean(distvals[distvals<Inf]),median(distvals))
  }
  overlapping[,(c("distmean","distmedian")):=rbindlist(pblapply(seq,calc,cl=threads))]
  result <- unnest_dt(data.table(row=seq_along(entries),nam=entries),"nam")
  result <- merge(result, overlapping)
  result[,list(list(nam[[which.min(distmedian)]])),keyby=row]$V1
}

#save.image()

for (markername in selmarker){
  marker_selname <- paste0(markername, "_selected")
  proteomes[, c(marker_selname):=get(markername)]
  proteomes[sapply(get(marker_selname),length)>1, c(marker_selname):=selectCandidate(get(marker_selname), alignmentfiles[[markername]],threads)]
  proteomes[, c(marker_selname):=sapply(get(marker_selname),function(el)ifelse(length(el)==0,NA,unlist(el)))]
}

#save.image()

getSequences <- function(markername, ids){
  al <- read.alignment(alignmentfiles[[markername]],"fasta")
  al <- as.data.table(unclass(al)[c("nam","seq")], key="nam")
  lengths <- unique(sapply(al$seq,nchar))
  if(length(lengths)>1)
    stop("Different lengths in alignment!?")
  ids[is.na(ids)] <- "<NOPE>"
  al <- rbind(al, data.table(nam="<NOPE>",seq=paste0(rep("-",lengths),collapse = "")))
  setkey(al,"nam")
  al[ids,seq]
}

proteomes[included==T,seq:=do.call(paste0,pblapply(selmarker,function(markername)getSequences(markername,get(paste0(markername, "_selected"))),cl=threads))]

#save.image()

fwrite(marker_stats,snakemake@output$stats,sep="\t")
proteomes[,(names(which(sapply(proteomes,is.list)))):=lapply(.SD, function(col)lapply(col,function(el)paste0(el,collapse = ","))),.SD=names(which(sapply(proteomes,is.list)))]
fwrite(proteomes,snakemake@output$proteomes,sep="\t")
fwrite(cbind(freq = 1, All = "All", proteomes[included==T, taxcols, with = F]),snakemake@output$krona_incl,sep="\t")
fwrite(cbind(freq = 1, All = "All", proteomes[included==F, taxcols, with = F]),snakemake@output$krona_excl,sep="\t")
