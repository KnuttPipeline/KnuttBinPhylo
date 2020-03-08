options(warn=2)

suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pbapply))


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


markerfile <- snakemake@input$markers
proteomes_file <- snakemake@input$proteomes
markersaligned <- snakemake@input$markeraligned
outputfile <- snakemake@output$tabular
outputfile_fas <- snakemake@output$seq
partfile <- snakemake@output$part
load(snakemake@input$tax)
threads <- snakemake@threads
sep <- "..."
refgroupsize <- 500
closerefs <- 5
distantrefs <- 15
marker_compl <- snakemake@params$markercompl
names(markersaligned) <- snakemake@params$markers
setDTthreads(1)
 
marker_stats <- fread(markerfile)
selmarker <- marker_stats[included==T, as.character(marker)]
proteomes <- fread(proteomes_file,select=c("Proteome","included","Organism ID",selmarker,paste0(selmarker,"_selected")))

readBinEntries <- function(markername){
  al <- read.alignment(markersaligned[[markername]],"fasta")
  al <- as.data.table(unclass(al)[c("nam","seq")],key="nam")
  al <- al[!unlist(strsplit(proteomes[[markername]],",",fixed=T))]
  al[,(c("bin","marker","id")):=tstrsplit(as.character(nam),sep,fixed=T)]
  al
}

binentries <- rbindlist(pblapply(selmarker, readBinEntries, cl=threads))
bins <- dcast(binentries, "bin~marker", value.var="nam", fun.aggregate=list, fill=list())
bins[,completion_perc:=100*rowSums(as.data.frame(lapply(.SD, function(col)sapply(col,function(el)length(el)>0))))/length(.SD),.SD=selmarker]
bins[,multicopy_perc:=100*rowSums(as.data.frame(lapply(.SD, function(col)sapply(col,length)>1)))/length(.SD),.SD=selmarker]
bins[,included:=completion_perc>marker_compl]

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
  result <- merge(result, overlapping,by="nam")
  result[,list(list(nam[[which.min(distmedian)]])),keyby=row]$V1
}

for (markername in selmarker){
  marker_selname <- paste0(markername, "_selected")
  bins[, c(marker_selname):=get(markername)]
  bins[sapply(get(marker_selname),length)>1, c(marker_selname):=selectCandidate(get(marker_selname), markersaligned[[markername]],threads)]
  bins[, c(marker_selname):=sapply(get(marker_selname),function(el)ifelse(length(el)==0,NA,unlist(el)))]
}

getSequences <- function(markername, ids){
  al <- read.alignment(markersaligned[[markername]],"fasta")
  al <- as.data.table(unclass(al)[c("nam","seq")], key="nam")
  lengths <- unique(sapply(al$seq,nchar))
  if(length(lengths)>1)
    stop("Different lengths in alignment!?")
  ids[is.na(ids)|ids==""] <- "<NOPE>"
  al <- rbind(al, data.table(nam="<NOPE>",seq=paste0(rep("-",lengths),collapse = "")))
  setkey(al,"nam")
  al[ids,seq]
}

proteomes <- proteomes[included==T,]
bins <- bins[included==T,]
setkey(bins,"bin")

proteomes[,seq:=do.call(paste0,pblapply(selmarker,function(markername)getSequences(markername,get(paste0(markername, "_selected"))),cl=threads))]
proteomes <- proteomes[,.(Proteome,`Organism ID`,seq)]
setkey(proteomes,"Proteome")
bins[,seq:=do.call(paste0,pblapply(selmarker,function(markername)getSequences(markername,get(paste0(markername, "_selected"))),cl=threads))]
bins[,seq:=pblapply(seq,s2c,cl=threads)]

bins[, distvals:=lapply(seq,function(binseq)unlist(pblapply(proteomes[,seq],function(refseq)JC(binseq,s2c(refseq)),cl=threads)))]
bins[, refgroup:=lapply(distvals,function(distvalsbin)proteomes$Proteome[distvalsbin %in% sort(distvalsbin)[1:refgroupsize]])]
bins[, refgroupdistvals:=lapply(distvals,function(distvalsbin)distvalsbin[distvalsbin %in% sort(distvalsbin)[1:refgroupsize]])]
bins[, refgroupdistvalorder:=lapply(refgroupdistvals,function(distvalsbin)order(distvalsbin))]
bins[, refgroupdistvals:=mapply(function(vals,order)vals[order],refgroupdistvals,refgroupdistvalorder,SIMPLIFY=F)]
bins[, refgroup:=mapply(function(vals,order)vals[order],refgroup,refgroupdistvalorder,SIMPLIFY=F)]
bins[, closerefgroup:=lapply(refgroup,function(group)group[1:closerefs])]
bins[, taxid:=sapply(closerefgroup,function(group)mrcatax(proteomes[unlist(group),`Organism ID`]))]
tax <- lookup(bins$taxid,threads=threads)
tax[,query:=NULL]
bins <- cbind(bins,tax)

#x=merge(unnest_dt(bins[,.(bin,closerefgroup)],closerefgroup),proteomes[,.(Proteome,`Organism ID`)],by.x="closerefgroup",by.y="Proteome")
#x=cbind(x,lookup(x$`Organism ID`))
#setorder(x,"bin")
#fwrite(x,"~/Schreibtisch/refs.tsv",sep="\t")

selectFarReferences <- function(binseq,refgroup,refgroupdistvals){
  refgroup <- unlist(refgroup)
  refgroupdistvals <- unlist(refgroupdistvals)
  seqs <- c(pblapply(proteomes[refgroup,seq],s2c,cl=threads))
  names(seqs) <- refgroup
  names(refgroupdistvals) <- refgroup
  result <- as.data.table(t(combn(names(seqs),2)))
  result$dist <- pbapply(as.matrix(result),1,function(i)JC(seqs[[i[[1]]]],seqs[[i[[2]]]]),cl=threads)
  result_mat <- matrix(NA, nrow=length(seqs), ncol=length(seqs), dimnames=list(names(seqs), names(seqs)))
  result_mat[as.matrix(result[, 1:2])] <- result$dist
  result_mat[lower.tri(result_mat)] <- t(result_mat)[lower.tri(result_mat)]
  result_mat <- as.dist(result_mat)
  clust <- hclust(result_mat,method="ward.D2")
  clusters <- cutree(clust,k=distantrefs)
  clusters <- data.table(proteome=names(clusters),group=clusters)
  clusters$disttobin <- refgroupdistvals[clusters$proteome]
  selected <- clusters[,.(proteome[[which.min(disttobin)]],min(disttobin)),by=group][,2:3]
  list(selected[[1]])
}
bins[, distantrefgroup:=mapply(selectFarReferences,seq,refgroup,refgroupdistvals)]

bins[,(names(which(sapply(bins,is.list)))):=lapply(.SD, function(col)lapply(col,function(el)paste0(el,collapse = ","))),.SD=names(which(sapply(bins,is.list)))]
fwrite(bins,outputfile,sep="\t")

selectedproteomes <- unique(unlist(c(strsplit(unlist(bins$distantrefgroup),",",fixed=T),strsplit(unlist(bins$closerefgroup),",",fixed=T))))
seqs <- proteomes[selectedproteomes,seq]
seqs <- append(seqs,gsub(",","",bins$seq,fixed=T))
nams <- c(selectedproteomes,bins$bin)

write.fasta(lapply(seqs,s2c),nams,outputfile_fas)


getLength <- function(markername)nchar(unclass(read.alignment(markersaligned[[markername]],"fasta"))$seq[[1]])
pos <- data.table(marker=selmarker, len=pbsapply(selmarker,getLength,cl=threads))
pos[,end:=cumsum(len)]
pos[,start:=shift(end,fill=0)+1]
pos[,model:=paste0("LG, ",marker,"=",as.character(start),"-",as.character(end))]
write(pos$model,partfile,ncolumns = 1)
