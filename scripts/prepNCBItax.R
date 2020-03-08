#!/usr/bin/env Rscript

##
## prepNCBItax.R - Create helper data and functions from the NCBI tax
##
## Knutt.org/Knutt2Reads2Bins

# From the names and nodes file this script produces the following objects
# in an RData file:
#
# tax_graph:
# Igraph object of the taxids as nodes and parent/child relationships
# as the edges. Directed from child to parent.
# tax_names:
# The taxid and name (scientific name) and rank of every 
# entry in a data.table
# lookup:
# See source code below
# mrcatax:
# See source code below
#

options(warn=2)

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(parallel))

name_file <- "reference_data/ncbitax/names.dmp"
node_file <- "reference_data/ncbitax/nodes.dmp"
out_file <- "test.RData"
threads <- detectCores()

if (exists("snakemake")) {
  name_file <- snakemake@input[["names"]]
  node_file <- snakemake@input[["nodes"]]
  out_file <- snakemake@output[[1]]
  threads <- snakemake@threads
}


options(mc.cores = threads)
setDTthreads(threads)


tax_names <- fread(
  file = name_file,
  sep = "|",
  header = F,
  col.names = c("taxid", "name", "unqiuename", "class", "wut")
)
# There is a weird column
tax_names[, wut := NULL]
cols <- c("name", "unqiuename", "class")
tax_names[, (cols) := mclapply(.SD, function(x)
  gsub('^\\s+|\\s+$', '', x)), .SD = cols]
tax_names <- tax_names[class == "scientific name", .(taxid, name)]


tax_nodes <-
  fread(
    node_file,
    sep = "|",
    header = F,
    col.names = c(
      "taxid",
      "parenttaxid",
      "rank",
      "embl",
      "division",
      "divflag",
      "geneticcode",
      "gc",
      "mitochgeneticcode",
      "mgcflag",
      "hidden",
      "hidden_subtree",
      "comments",
      "wut"
    )
  )
tax_nodes[, wut := NULL]
cols <- c("rank", "embl", "comments")
tax_nodes[, (cols) := mclapply(.SD, function(x)
  gsub('^\\s+|\\s+$', '', x)), .SD = cols]
tax_graph = graph_from_edgelist(as.matrix(tax_nodes[, c(1, 2)]))


setkey(tax_names, taxid)
setkey(tax_nodes, taxid)
tax_names <- tax_names[tax_nodes]
tax_names <- tax_names[, .(name, rank), keyby = taxid]
rm(tax_nodes)

# Returns the taxonomy as a data.table with query being the given taxids,
# each requested level as column. If filling is activated, gaps will be
# filled with "(Taxonomy Gap)" and missing information with "Unclassified"
# Multithreading is done with pbapply
#
# Needs tax_graph and tax_names
lookup <-
  function(taxids,
           taxlevels = c("superkingdom",
                         "phylum",
                         "class",
                         "order",
                         "family",
                         "genus",
                         "species"),
           threads = 1,
           fillgaps = T) {
    suppressPackageStartupMessages(require(igraph))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(pbapply))
    uniquetaxids <- unique(taxids)
    df <- rbindlist(pblapply(uniquetaxids, function(x) {
      if (x %in% tax_names$taxid && x != 1)
        data.table(query = x, taxid = as.numeric(all_simple_paths(tax_graph, x, 1)[[1]]))
      else
        data.table(query = x, taxid = x)
    }, cl = threads))
    setkey(df, taxid)
    df <- tax_names[df]
    df <-
      data.table::dcast(df[rank != "no rank", ], query ~ rank, value.var = "name")
    if (length(taxlevels) != 0) {
      missingcols <- taxlevels[!taxlevels %in% colnames(df)]
      if (length(missingcols) > 0)
        df[, c(missingcols) := NA]
      df <- df[, c("query", taxlevels), with = FALSE]
    }
    
    tax <- df[data.table(query = taxids), on = "query"]
    if (fillgaps) {
      for (taxi in rev(seq_along(taxlevels)[-1])) {
        currentlevel <- taxlevels[taxi]
        upperlevel <- taxlevels[taxi - 1]
        tax[!is.na(get(currentlevel)) &
              is.na(get(upperlevel)), (upperlevel) := "(Taxonomy Gap)"]
      }
      tax[is.na(tax)] <- "Unclassified"
    }
    tax
  }

# Returns the first (children -> parents) common ancestor for the
# given taxids. 1 if there isn't one (root)
mrcatax <- function(taxids) {
  suppressPackageStartupMessages(require(igraph))
  suppressPackageStartupMessages(require(data.table))
  uniquetaxids <- unique(taxids)
  if (length(uniquetaxids) == 0)
    return(NA)
  lineages = lapply(uniquetaxids, function(x)
    as.data.frame(matrix(rev(
      as.numeric(all_simple_paths(tax_graph, x, 1)[[1]])
    ), nrow = 1)))
  lineages <- rbindlist(lineages, fill = T)
  lineage_overlapping <- which(apply(lineages, 2, function(col)
    length(unique(col)) == 1))
  lineage_overlapping <- tail(lineage_overlapping, 1)
  lineages[1, lineage_overlapping, with = F][[1]]
}

save(tax_graph, tax_names, lookup, mrcatax, file = out_file)
