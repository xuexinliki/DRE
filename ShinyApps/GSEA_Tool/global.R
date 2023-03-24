suppressPackageStartupMessages(library(BiocManager))
options(repos = BiocManager::repositories())
# suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyscreenshot))
suppressPackageStartupMessages(library(ggplot2))
# devtools::install_github("nicolash2/gggsea")
# suppressPackageStartupMessages(library(gggsea))
# suppressPackageStartupMessages(library(ggrepel))
# suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(shinydlplot))
suppressPackageStartupMessages(library(shinyjs))
# suppressPackageStartupMessages(library(randomcoloR))
suppressPackageStartupMessages(library(fgsea))
# suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(dplyr))

# helper functions --------------------------------------------------------

# data --------------------------------------------------------------------
organism <- readRDS("DB/Organism_Names.RDS")
organism <- list("Main Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T)],
                 "Other Species" = organism[grep("Rattus|Musculus|Homo|Dresophila|Danio", organism, ignore.case = T, invert = T)])
cclasses <- c("H", paste("C",1:8, sep = ""))

# registerDoSEQ()



