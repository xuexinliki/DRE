library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(limma)
library(data.table)

files <- list.files(path = "MSigDB/", pattern = ".RDS", ignore.case = T, full.names = T)
files <- sort(files)

cmap <- readRDS("DREMIT/data/D1_short_t_matrix_12434_4690_LINCS.rds")
row.names(cmap) <- toupper(row.names(cmap))
ref <- readRDS("DREMIT/data/sig_id_table_LINCS_short.rds")
corganism <- gsub(".*MSigDB_Species_(.*).RDS","\\1",files[j], ignore.case = T)
print(corganism)
cdir <- paste("GSEA_Result/",gsub("\\s+","_",corganism),"/",sep = "")
dir.create(cdir, recursive = T)
for(i in 1:ncol(cmap)){
crank <- NULL
crank <- cmap[,i]
names(crank) <- row.names(cmap)
crank = sort(crank, decreasing = TRUE)

for(j in 1:length(files)){
sigsets <- NULL
sigsets <- readRDS(files[j])
sigsets <- lapply(sigsets, function(x){x <- toupper(x[!is.na(x)])})
cresult <- NULL
cresult <- fgseaMultilevel(sigsets, stats = crank, nPermSimple = 10000, eps = 0)
cresult <- cresult[order(cresult$padj, decreasing = F),]
cresult <- cresult[which(cresult$padj < 0.05),]
cresult <- data.frame(Signature = colnames(cmap)[i],
Organism = gsub("(.*?)\\|(.*?)\\|(.*)$","\\1",cresult$pathway, ignore.case = T),
Category = gsub("(.*?)\\|(.*?)\\|(.*)$","\\2",cresult$pathway, ignore.case = T),
cresult)
saveRDS(cresult, paste(cdir, corganism, "_", colnames(cmap)[i], ".RDS", sep = ""))
}
}

print("Completed!")
