library(reshape2)

files <- list.files(path = "GSEA_Result/")
files <- sort(files)
for(i in 1:length(files)){
print(i)
cname <- NULL
cname <- gsub("_GSEA_PAdj005.RDS","",files[i], ignore.case = T)
x <- NULL
x <- readRDS(files[i])
corganism <- NULL
corganism <- unique(gsub("^(.*?)\\|.*","\\1",x$pathway[1], ignore.case = T))

data <- NULL
data <- dcast(x, pathway~Signature, value.var = "ES")
row.names(data) <- data$pathway
data <- data[,which(colnames(data) != "pathway")]

# The following codes are adapted from DREIMT (http://www.dreimt.org/) with partial modification to fit our data
taufun <- function(NormScore, DrugCol){
  x <- sum(abs(DrugCol) < abs(NormScore))/length(DrugCol) * 100
  return(x)
}
data.pos <- NULL
data.neg <- NULL
data.pos <- data
data.neg <- data

data.pos[data.pos <= 0] <- NA
data.neg[data.neg > 0] <- NA

num <- NULL
not.na <- NULL
den <- NULL
num <- outer(rowSums(data.pos, na.rm = TRUE), colSums(data.pos, na.rm = TRUE), "+") - data.pos
not.na <- !is.na(data.pos)
den <- outer(rowSums(not.na), colSums(not.na), "+") - 1
norm.fact <- NULL
norm.pos <- NULL
norm.fact <- num/den
norm.pos <- data.pos/norm.fact

num <- NULL
not.na <- NULL
den <- NULL
num <- outer(rowSums(data.neg, na.rm = TRUE), colSums(data.neg, na.rm = TRUE), "+") - data.neg
not.na <- !is.na(data.neg)
den <- outer(rowSums(not.na), colSums(not.na), "+") - 1
norm.fact <- NULL
norm.neg <- NULL
norm.fact <- num/den
norm.neg <- data.neg/abs(norm.fact)

temp <- NULL
temp <- melt(data.frame(pathway = row.names(norm.neg), norm.neg))
summary(temp$value[!is.na(temp$value)])

norm.ES <- NULL
norm.ES <- melt(data.frame(pathway = row.names(norm.pos), norm.pos))
temp <- NULL
temp <- melt(data.frame(pathway = row.names(norm.neg), norm.neg))
norm.ES <- rbind(norm.ES[which(!is.na(norm.ES$value)),], temp[which(!is.na(temp$value)),])
norm.ES <- data.frame(Organism = corganism, norm.ES)
colnames(norm.ES) <- c("Organism","pathway","Signature","Norm_ES")
saveRDS(norm.ES, paste("GSEA_Norm_ES_Padj005/", cname, "_Norm_ES_Padj005.RDS", sep = ""))
  
  norm.ES <- dcast(norm.ES, pathway ~ Signature, value.var = "Norm_ES")
  print(dim(norm.ES))
  row.names(norm.ES) <- norm.ES$pathway
  norm.ES <- norm.ES[,which(colnames(norm.ES) != "pathway")]
  norm.ES[is.na(norm.ES)] <- 0
  ES.tau <- NULL
  for(j in 1:ncol(norm.ES)){
    print(j)
    current <- NULL
    cref <- NULL
    cref <- norm.ES[which(norm.ES[,j] != 0),j]
    names(cref) <- row.names(norm.ES[which(norm.ES[,j] != 0),])
    current <- sapply(cref, taufun, cref)
    current <- data.frame(Signature = colnames(norm.ES)[j], pathway = names(cref), Tau = current*sign(cref))
    ES.tau <- rbind(ES.tau, current)
  }
  
  saveRDS(ES.tau, paste("GSEA_Tau_Padj005/", cname, "_GSEA_Tau_Padj005.RDS", sep = ""))

  x$Tau_Score <- ES.tau[match(paste(x$pathway,x$Signature), paste(ES.tau$pathway,ES.tau$Signature)),"Tau"]
  saveRDS(x, paste("../GSEA_Result_Padj005_Tau/",gsub("GSEA_PAdj005.RDS","GSEA_PAdj005_Tau.RDS",files[i], ignore.case = T), sep = ""))
  
}