suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(igraph))


colors <- tolower(c("#E66543","#4133E9","#CD79D5","#4AA053","#071013"))

tree_file <- "output/tree/RAxML_bipartitionsBranchLabels.run"
bin_file <- "output/tree/totree.tsv"
reference_file <- "reference_data/proteomes.tsv"
tax_file <- "reference_data/ncbi_tax/ncbi_tax.RData"

leaf_type_file <- "leaf_types.txt"
label_file <- "species_labels.txt"
bold_label_file <- "bold_labels.txt"

tree <- read.tree(tree_file)
bins <- fread(bin_file)

types <- data.table(id=tree$tip.label, key="id")
types[,isclose:=-1]
types[,isdistant:=-1]
closerefs <- unique(unlist(strsplit(bins$closerefgroup, ",", fixed=T)))
distantrefs <- unique(unlist(strsplit(bins$distantrefgroup, ",", fixed=T)))
types[!id %in% bins$bin, isclose:=as.numeric(id %in% closerefs)]
types[!id %in% bins$bin, isdistant:=as.numeric(id %in% distantrefs)]

write(paste0("DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tLeaf Type\nCOLOR\t",colors[[1]],"\nFIELD_SHAPES\t1\t1\nFIELD_LABELS\tClose reference\tDistant reference\nSHOW_LABELS\t1\nFIELD_COLORS\t",colors[[1]],"\t",colors[[2]],"\nDATA\n"),leaf_type_file)
fwrite(types, leaf_type_file, append=T, sep="\t")

tree_igraph <- as.igraph.phylo(tree,T)
load(tax_file)
write("LABELS\nSEPARATOR TAB\nDATA\n",label_file)
refs <- fread(reference_file,key="Proteome",select=c("Proteome", "Organism ID", "BUSCO"))
refs <- refs[unique(c(closerefs,distantrefs))]
r <- "([0-9]*\\.[0-9]+|[0-9]+)%"
regex <- paste0("C:",r,"\\[S:",r,",D:",r,"\\],F:",r,",M:",r,",n:(\\d+)")
refs[BUSCO!="",BUSCO:=sub(regex,"\\1,\\2,\\3,\\4,\\5,\\6",BUSCO)]
refs[BUSCO!="",c("CompletionTotal","CompletionSingle","CompletionDuplicated","Fragmented","Missing","ProteinCount"):=tstrsplit(BUSCO,",",T)]
setkey(refs,"Proteome")
mrca <- function(taxids){
    taxids <- unique(taxids[!is.na(taxids)])
    if(length(taxids)==0){
        return(NA)
    }
    if(length(taxids)==1){
        return(taxids)
    }
    lineages <- sapply(taxids,function(tax)rev(as.numeric(all_simple_paths(tax_graph,tax,1)[[1]])))
    shortest <- min(sapply(lineages,length))
    taxid <- NA
    for(step in seq(shortest)){
        currentlevels <- unique(sapply(lineages, function(lin)lin[[step]]))
        if(length(currentlevels)>1)
            return(taxid)
        taxid <- currentlevels
    }
    return(taxid)
}
rootnodes <- names((which(degree(tree_igraph, mode="in")==0)))
for(rootnode in rootnodes){
    for(node in rev(names(bfs(tree_igraph,rootnode)$order))){
        if(node %in% refs$Proteome){
            taxid <- refs[node, `Organism ID`]
        }else if (node %in% bins$bin) {
           taxid <- NA
        }else{
           subnodes <- neighbors(tree_igraph, v=node, "out")
           subids <- vertex_attr(tree_igraph, "taxid", subnodes)
           taxid <- mrca(subids)
        }
        tree_igraph <- set_vertex_attr(tree_igraph, "taxid", node, taxid)
    }
}
getLCAID <- function(node){
    if(degree(tree_igraph, node, "out")==0)
        return(node)
    children <- names(neighbors(tree_igraph,node,"out"))
    children_leafoffspring <- lapply(children, function(child)subcomponent(tree_igraph, child, "out"))
    children_leafoffspring <- sapply(children_leafoffspring, function(childoff)names(childoff[degree(tree_igraph, childoff, "out")==0][[1]]))
    return(paste0(children_leafoffspring[[1]],"---",children_leafoffspring[[2]]))
}
labels <- data.table(node=names(V(tree_igraph)),label=sapply(names(V(tree_igraph)),getLCAID), taxid=vertex_attr(tree_igraph, "taxid"))
labels <- merge(labels,tax_names,by="taxid")
labels[label %in% refs$Proteome,name:=paste0(name," (",label,")")]
fwrite(labels[,.(label,name)], label_file, append=T, sep="\t")

write(paste0("DATASET_STYLE\nSEPARATOR TAB\nDATASET_LABEL\tBins\nCOLOR\t",colors[[2]],"\nDATA\n"),bold_label_file)
bolds <- data.table(id=tree$tip.label, key="id")
bolds <- bolds[id %in% bins$bin, ]
bolds[,(c("what","color","size","style","bgcolor")):=list(list("label"),list("node"),list("#000000"),list("1"),list("bold-italic"))]
fwrite(bolds, bold_label_file, append=T, sep="\t")

