#/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(dplyr))

totreefile <- snakemake@input$totree
treefile <- snakemake@input$tree
protfile <- snakemake@input$proteomes
gtdbfile <- snakemake@input$gtdb
pruningdir <- snakemake@output$prunedir
dir.create(pruningdir)
taxdir <- snakemake@output$taxdir
dir.create(taxdir)
treecolors <- snakemake@output$colors
labels <- snakemake@output$labels
popups <- snakemake@output$popups
outtreefile <- snakemake@output$tree
markers <- snakemake@params$markers

taxcols <- c("superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species")


totree <- fread(totreefile)
totree[,distantrefgroup:=strsplit(distantrefgroup,",",T)]
totree[,closerefgroup:=strsplit(closerefgroup,",",T)]
distantrefs <- unique(unlist(totree$distantrefgroup))
closerefs <- unique(unlist(totree$closerefgroup))
setkey(totree,bin)

treeall <- read.raxml(treefile)
tree <- treeall@phylo
tree$node.label <- as.data.table(as_tibble(treeall))[,.(parent,node,bootstrap)][setnames(as.data.table(tree$edge),c("parent","node")),on=c("parent","node")]$bootstrap

proteomes <- fread(protfile,select=c("Proteome", paste0(markers,"_selected"),taxcols,"Genome assembly ID"))
setnames(proteomes, "Genome assembly ID", "assembly")
setkey(proteomes, Proteome)
proteomes <- proteomes[tree$tip.label,nomatch=0]
gtdb <- fread(gtdbfile, header=F, col.names=c("assembly","GTDB_tax"),sep="\t")
gtdb[,assembly:=sub("^.._(.+)$","\\1",assembly)]
gtdb[,GTDB_tax:=gsub(";","; ",GTDB_tax)]
proteomes <- gtdb[proteomes, on="assembly"]
rm(gtdb)
setkey(proteomes, Proteome)
proteomes[,close:=Proteome %in% closerefs]
proteomes[,distant:=Proteome %in% distantrefs]
tree_edgedf <- cbind(data.table(From=tree$edge[,1], To=tree$edge[,2]), edgelen=tree$edge.length)
setorder(tree_edgedf,edgelen)

getLCAID <- function(tree, nodeid){
    lcanodeids <- sapply(Descendants(tree, Children(tree, nodeid), type="tips"),"[[",1)
    if(length(lcanodeids)==0){
        lcanodeids <- nodeid
    }
    paste0(tree$tip.label[lcanodeids],collapse="|")
}
getTaxAssigment <- function(tree, nodeid, proteomes){
    tips <- tree$tip.label[Descendants(tree, nodeid, type="tips")[[1]]]
    nodetax <- proteomes[tips, apply(.SD, 2, function(col)list(unique(col))),nomatch=0,.SD=taxcols]
    if(length(tips)==0)
        return(data.table(tipcount=nrow(nodetax)))
    taxcount <- sapply(nodetax,function(x)length(unlist(x)))
    selectedlevel <- last(names(which(taxcount==min(taxcount))))
    cbind(data.table(tipcount=length(tips),level=selectedlevel, taxvals=paste0(unlist(nodetax[[selectedlevel]]),collapse=", "), taxcount=length(unlist(nodetax[[selectedlevel]]))),nodetax[,lapply(.SD,unlist),,.SD=names(which(taxcount==1))])
}

#outgroup_node <- tree_edgedf[edgelen==max(edgelen),To]
#outgroup_nodes <- Descendants(tree, outgroup_node, type="tips")[[1]]
#ingroup_nodes <- setdiff(1:tree$Nnode, outgroup_nodes)
#tree <- root(tree, outgroup_nodes)
write.tree(tree, outtreefile)
taxlabels <- rbindlist(lapply(unique(unlist(as.list(tree$edge))),function(i)cbind(node=i,getTaxAssigment(tree=tree,i,proteomes=proteomes))), fill=T)
setkey(taxlabels,node)
taxlabels[tipcount>1,id:=getLCAID(tree,first(node)),by=node]



headerlines <- c("TREE_COLORS","SEPARATOR TAB","DATA")
#outgroupline <- paste0(c(getLCAID(tree,outgroup_node),"range","#b3c6ff","Automatic outgroup"),collapse="\t")
binlines <- paste0(totree[bin %in% tree$tip.label, bin],"\tlabel\t#000000\tbold")
proteomes[close==T&distant==F,color:="#fb8b24"]
proteomes[close==T&distant==T,color:="#5f0f40"]
proteomes[close==F&distant==T,color:="#0f4c5c"]
proteomelines <- proteomes[,paste0(Proteome, "\tlabel\t",color)]
fileConn<-file(treecolors)
writeLines(c(headerlines,binlines,proteomelines), fileConn)
close(fileConn)

headerlines = c("POPUP_INFO","SEPARATOR TAB","DATA")
popupsdat <- proteomes[,paste0(Proteome,"\t","Proteome Info","\t","<h1>",Proteome,"</h1>",
    "<p>NCBI taxonomy: <b>", apply(.SD,1,paste,collapse="; "),"</b></p>",
    "<p>GTDB taxonomy: <b>", GTDB_tax,"</b></p>",
    "<p>UniProt: <b><a target='_blank' href='https://www.uniprot.org/proteomes/",Proteome,"'>",Proteome,"</a></b></p>",
    "<p>Assembly: <b><a target='_blank' href='https://www.ncbi.nlm.nih.gov/assembly/",assembly,"'>",assembly,"</a></b></p>")
 , .SD=taxcols]
fileConn<-file(popups)
writeLines(c(headerlines,popupsdat), fileConn)
close(fileConn)

headerlines <- c("LABELS","SEPARATOR TAB","DATA")
proteomelabels <- proteomes[,paste0(Proteome,"\t",Proteome," (",species,")")]
nodelabels <- taxlabels[!is.na(id),paste0(id,"\t",taxvals)]
binlabels <- paste0(totree[bin %in% tree$tip.label, bin],"\t", totree[bin %in% tree$tip.label, bin])
fileConn<-file(labels)
writeLines(c(headerlines,proteomelabels, nodelabels, binlabels), fileConn)
close(fileConn)

writePruneFile <- function(bin){
    fileConn <- file(file.path(pruningdir,paste0(bin,"_prune.txt")))
    writeLines(c("PRUNE", "DATA", bin,unlist(totree[bin,unique(unlist(c(distantrefgroup, closerefgroup)))])), fileConn)
    close(fileConn)
}
for(bin in unique(totree$bin)) writePruneFile(bin)

color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
writeTaxFile <- function(leveltouse){
    curtaxlabels <- taxlabels[!is.na(get(leveltouse)),node,by=leveltouse]
    curtaxlabels[,nodes:=list(list(node)),by=leveltouse]
    curtaxlabels[,parent:=Ancestors(tree, node, type = "parent")]
    tomark <- curtaxlabels[apply(curtaxlabels,1,function(row)!row$parent %in% row$nodes),.(get(leveltouse),node)]
    tomark <- tomark[,.(lca=getLCAID(tree,first(node)), V1),by=node]
    colorsel <- tolower(col2hex(sample(color, length(unique(curtaxlabels[[leveltouse]])))))
    tomark[,color:=colorsel[[.GRP]],by=V1]
    rows <- tomark[,paste0(lca,"\trange\t",color,"\t",V1)]
    fileConn <- file(file.path(taxdir,paste0(leveltouse,"_marked.txt")))
    writeLines(c("TREE_COLORS", "SEPARATOR TAB", "DATA", rows), fileConn)
    close(fileConn)
  #  rows <- tomark[,paste0(lca,"\t",V1,"\t-1\t",color,"\tnormal\t2")]
  #  fileConn <- file(file.path(taxdir,paste0(leveltouse,"_marked_text.txt")))
  #  writeLines(c("DATASET_TEXT", "SEPARATOR TAB", paste0("DATASET_LABEL\t",leveltouse),"COLOR\t#ff0000", "DATA", rows), fileConn)
  #  close(fileConn)
}
for(tax in taxcols) writeTaxFile(tax)
