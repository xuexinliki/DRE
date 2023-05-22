suppressPackageStartupMessages(library(BiocManager))
options(repos = c(BiocManager::repositories(),"https://bioconductor.org/packages/release/bioc/src/contrib/"))
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyscreenshot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(shinydlplot))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))

# helper functions --------------------------------------------------------

# data --------------------------------------------------------------------
ref <- readRDS("DB/sig_id_table_LINCS_short.rds")
cmap <- readRDS(url("https://www.dropbox.com/s/3uacuq4n52lf1ia/Filtered_D1_short_t_matrix_12434_4690_LINCS.rds?dl=1","rb"))
row.names(cmap) <- toupper(row.names(cmap))
compoundurl <- "https://pubchem.ncbi.nlm.nih.gov/compound/"
organism <- readRDS("DB/Organism_Names.RDS")
organism <- list("Main Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T)],
                 "Other Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T, invert = T)])
classcols <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold","#a65628","#999999","#9e2a2b","grey","#58bc82","purple")
groupcols <- c("#5F75A5","#91D376","#D75B58","#F5BFD3","#A8C9EA","#B09C16","#F69F99","#AC79A3","#E89314","#EAD256","#78706E","#D1A5CA","#F7C277","#569794","#B9B0AC","#99785E","#5FA346","#8DBCB6","#CC7296","#D3B6A5")
fillcols <- hcl.colors(10, palette = "Spectral")
reflinks <- readRDS("DB/DB_Links.RDS")



