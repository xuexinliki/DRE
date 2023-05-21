suppressPackageStartupMessages(library(BiocManager))
options(repos = BiocManager::repositories())
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyscreenshot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(shinydlplot))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(randomcoloR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(qusage))
suppressPackageStartupMessages(library(GeneOverlap))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))

# helper functions --------------------------------------------------------

# data --------------------------------------------------------------------
organism <- readRDS("DB/Organism_Names.RDS")
organism <- list("Main Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T)],
                 "Other Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T, invert = T)])
orgurl <- readRDS("DB/Human_Gene_Symbol_URL.RDS")

cclasses <- c("H", paste("C",1:8, sep = ""))
fillcols <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))

# registerDoSEQ()



