suppressPackageStartupMessages(library(BiocManager))
options(repos = c(BiocManager::repositories(),"https://bioconductor.org/packages/release/bioc/src/contrib/"))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinydlplot))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(ggparallel))
suppressPackageStartupMessages(library(parcoords))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggsankey))

organism <- readRDS("DB/Organism_Names.RDS")
organism <- list("Main Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T)],
                 "Other Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T, invert = T)])
drugname <- readRDS("DB/L1000_DB_Drug_Names.RDS")
classcols <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold","#a65628","#999999","#9e2a2b","grey","#58bc82","purple")
groupcols <- c("#5F75A5","#91D376","#D75B58","#F5BFD3","#A8C9EA","#B09C16","#F69F99","#AC79A3","#E89314","#EAD256","#78706E","#D1A5CA","#F7C277","#569794","#B9B0AC","#99785E","#5FA346","#8DBCB6","#CC7296","#D3B6A5")
taucols <- RColorBrewer::brewer.pal(11,"RdYlBu")
pathwayclass <- c("H", paste("C",1:8, sep = ""))
csummary <- readRDS("DB/MSigDB_URLs.RDS")
reflinks <- readRDS("DB/Human_Gene_Symbol_URL.RDS")
cmap <- readRDS(url("https://www.dropbox.com/s/z2hbx36fl6joryx/D1_short_t_matrix_12434_4690_LINCS.rds?dl=1","rb"))
# cmap <- readRDS("DB/D1_short_t_matrix_12434_4690_LINCS.rds")
ref <- readRDS("DB/sig_id_table_LINCS_short.rds")
ref$common_name <- stringr::str_to_title(ref$common_name)
ref$Drug_Name <- paste(ref$common_name, " [Profile ID=",ref$sig_id,"]", sep = "")
orgurl <- readRDS("DB/GSEA_Result_Padj005_Tau_Abs80_Mini_URL.RDS")

sca_heatmap <- function(plotx){
  if(nrow(plotx)>0){
    p <- plot_ly(x = colnames(plotx), y = rownames(plotx), z = as.matrix(plotx), type = "heatmap", colors = "RdYlBu", strokes = "black") %>%
      layout(xaxis = list(tickangle = -45))
    p
  }else{
    stop("Drug not found")
  }
}

# squish_trans from: https://stackoverflow.com/users/4303162/stibu
squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from) / factor
    x[ito] <- from + (to - from) / factor + (x[ito] - to)
    return(x)
  }
  
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from) / factor
    ito <- x >= from + (to - from) / factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from) / factor))
    return(x)
  }
  
  # return the transformation
  return(scales::trans_new("squished", trans, inv))
}

