suppressPackageStartupMessages(library(BiocManager))
options(repos = BiocManager::repositories())
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyscreenshot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(shinydlplot))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(dplyr))

# helper functions --------------------------------------------------------

# data --------------------------------------------------------------------
sigsets <- readRDS("DB/L1000_Dreimt_Consensus_Drug_Signatures.RDS")
sigsets <- lapply(sigsets, function(x){x <- toupper(x); return(x)})

