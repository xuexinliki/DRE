suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinydlplot))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(markdown))

organism <- readRDS("DB/Organism_Names.RDS")
organism <- list("Main Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T)],
                 "Other Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T, invert = T)])
classcols <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold","#a65628","#999999","#9e2a2b","grey","#58bc82","purple")
taucols <- RColorBrewer::brewer.pal(11,"RdYlBu")
csummary <- readRDS("DB/MSigDB_Disease_Subtypes_All_Species.RDS")
csummary <- csummary[which(csummary$SubCategory %in% c("CP:KEGG","CP:WIKIPATHWAYS")),]
diseases <- sort(unique(csummary$DiseaseType))
pathwayclass <- sort(unique(csummary$Category))

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
