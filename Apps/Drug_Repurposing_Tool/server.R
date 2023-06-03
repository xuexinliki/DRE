sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

function(input, output, session){

    goi <- reactive({
      req(input$corganism)
      req(input$run)
      req(input$goi)
      print("goi..")
        x <- unlist(strsplit(input$goi, split = ","))
        x <- unlist(strsplit(x, split = "\n"))
        x <- unlist(strsplit(x, split = " "))
        x <- x[which(!is.na(x))]
        x <- x[which(x != "")]
        return(x)
    })
    
    corganism <- reactive({
      input$corganism
    })
    
    data <- reactive({
      req(input$goi)
      req(goi())
      print("input$goi:")
      print(input$goi)
      
      showModal(modalDialog("Running analysis (whole analysis may take up to 15 minutes)..", footer=NULL))
        
        numCores<-detectCores()
        cl <- makeCluster(numCores/2)
        registerDoParallel(cl)

        print("GOI:")
        print(goi())
        goi <- goi()
        
        data <- NULL
        data <- foreach (i=1:ncol(cmap), .combine=rbind, .export = c("cmap")) %dopar% {
          print(i)
          crank <- NULL
          crank <- cmap[,i]
          names(crank) <- row.names(cmap)
          crank = sort(crank, decreasing = TRUE)
          cresult <- NULL
          print(goi)
          cresult <- fgsea::fgseaMultilevel(list(UserQuery=goi), stats = crank, nPermSimple = 200, eps = 0)
          cresult <- data.frame(cresult)
          cresult <- cresult[,which(!colnames(cresult) %in% c("leadingEdge","size"))]
          cresult <- cresult[order(cresult$padj, decreasing = F),]
          cresult <- data.frame(cresult)
          if(length(which(cresult$padj < 0.05)) > 1){
            cresult <- cresult[which(cresult$padj < 0.05),]
          }
          if(nrow(cresult) > 0){
            cresult <- data.frame(Signature = colnames(cmap)[i],
                                  cresult[,which(colnames(cresult) != "pathway")])
            return(cresult)
            # data <- rbind(data, cresult)
          }
        }
        
        print("Before stopCluster..")
        
        stopCluster(cl)
        registerDoSEQ()
        
        print("Done with GSEA..")
        
        data <- data.frame(data)
        if(length(which(data$padj < 0.05)) > 10){ # != nrow(data)
          data <- data[which(data$padj < 0.05),]
        }
        print("Data dim:")
        print(dim(data))
        print(head(data))
        if(nrow(data) > 2000){
          data <- data[order(data$pval, decreasing = F),]
          data <- data[1:ifelse(nrow(data) > 2000, 2000, nrow(data)),]
        }
        
        data$Drug_Name <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"common_name"])
        data$Drug_Name <- paste(stringr::str_to_title(data$Drug_Name)," [Profile ID=",gsub("sig_","",data$Signature, ignore.case = T),"]", sep = "")
        data$Mechanism <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"MOA"])
        data$Drug_Status <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"Status"])
        data$PubChem_ID <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"pubchem_cid"])
        data$LINCS_ID <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"lincs_id"])
        data[which(!is.na(data$PubChem_ID)),"PubChem_ID"] <- paste0("<a href='",compoundurl,data[which(!is.na(data$PubChem_ID)),"PubChem_ID"],"' target='_blank'>",data[which(!is.na(data$PubChem_ID)),"PubChem_ID"],"</a>")
        data[which(!is.na(data$LINCS_ID)),"LINCS_ID"] <- paste0("<a href='http://www.ilincs.org/ilincs/perturbagen/compound/",data[which(!is.na(data$LINCS_ID)),"LINCS_ID"],"' target='_blank'>",data[which(!is.na(data$LINCS_ID)),"LINCS_ID"],"</a>")
        data$Source_ID <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"source_id"])
        data$Drug_Specificity <- unlist(ref[match(gsub("sig_","",data$Signature, ignore.case = T), ref$sig_id),"DSS"])
        removeModal()

        showModal(modalDialog("Calculating tau score..", footer=NULL))
        
        taufun <- function(NormScore, DrugCol){
          x <- sum(abs(DrugCol) < abs(NormScore))/length(DrugCol) * 100
          return(x)
        }
        
        data$pathway <- "GOI"
        dcastdata <- NULL
        print("Data:")
        print(head(data))
        dcastdata <- dcast(data, pathway~Signature, value.var = "ES")
        row.names(dcastdata) <- dcastdata$pathway
        dcastdata <- dcastdata[,which(colnames(dcastdata) != "pathway")]
        dcastdata[is.na(dcastdata)] <- 0
        dcastdata <- as(as.matrix(dcastdata), "dgCMatrix")
        print("Post dgCMatrix..")
        clink <- NULL
        clink <- reflinks[match(paste(gsub(" ","_",corganism(), ignore.case = T), "_Norm_ES_DCast_Padj005.RDS", sep = ""), reflinks$File),"Link"]
        clink <- gsub("dl=0","dl=1",clink)
        dbES <- readRDS(url(clink,"rb"))
        dbES <- dbES[,which(colnames(dbES) %in% colnames(dcastdata))]
        dbES <- dbES[,match(colnames(dcastdata), colnames(dbES))]
        dcastdata <- rbind(dbES, dcastdata)
        rm(dbES)
        gc()
        
        data.pos <- NULL
        data.neg <- NULL
        data.pos <- dcastdata
        data.neg <- dcastdata
        rm(dcastdata)
        print("removed dcastdata..")
        
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
        rm(data.pos)
        
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
        rm(data.neg)
        rm(den)
        rm(not.na)
        rm(num)
        rm(norm.fact)
        
        norm.ES <- NULL
        norm.ES <- melt(data.frame(pathway = row.names(as.matrix(norm.pos)), as.matrix(norm.pos)))
        rm(norm.pos)
        temp <- NULL
        temp <- melt(data.frame(pathway = row.names(as.matrix(norm.neg)), as.matrix(norm.neg)))
        rm(norm.neg)
        norm.ES <- rbind(norm.ES[which(!is.na(norm.ES$value)),], temp[which(!is.na(temp$value)),])
        rm(temp)
        colnames(norm.ES) <- c("pathway","Signature","Norm_ES")
        gc()
        
        norm.ES <- dcast(norm.ES, pathway ~ Signature, value.var = "Norm_ES")
        norm.ES[is.na(norm.ES)] <- 0
        row.names(norm.ES) <- norm.ES$pathway
        norm.ES <- norm.ES[,which(colnames(norm.ES) != "pathway")]
        
        i <- which(row.names(norm.ES) == "GOI")
        print(row.names(norm.ES)[i])
        print("Calculating Tau..")
        
        for(j in 1:dim(norm.ES)[2]){
          norm.ES[i,j] <- ((sign(norm.ES[i,j]) * 100)/dim(norm.ES)[1]) * sum(abs(norm.ES[,j]) < abs(norm.ES[i,j]))
        }
        
        norm.ES <- norm.ES[which(row.names(norm.ES) == "GOI"),]
        gc()
        norm.ES <- melt(norm.ES)
        
        data$Tau <- norm.ES[match(data$Signature,norm.ES$variable),"value"]
        rm(norm.ES)
        if(length(which(abs(data$Tau) > 0.9)) > 1){
          data <- data[which(abs(data$Tau) > 0.9),]
        }
        print("After Tau:")
        print(head(data))
        data$FDR <- data$padj
        data <- data[,c("LINCS_ID","Drug_Name","Mechanism","Drug_Status","Tau","FDR","ES","Drug_Specificity","Source_ID","PubChem_ID")]
        data <- data.frame(data)
        data <- data[order(data$Tau, decreasing = F),]
        data <- data.frame(Effect = ifelse(data$Tau > 0, "Postive Tau", "Negative Tau"), data)
        data$Effect[which(data$Tau == 0)] <- "Zero Tau"
        data$Effect <- factor(data$Effect, levels = c("Postive Tau", "Negative Tau", "Zero Tau"))
        removeModal()
        return(data)
    })
    
    # output$contents <- renderPrint({
    #     data()
    # })
    
    PlotObject <- reactive({
        plotx <- data()
        plotx <- plotx[which(plotx$FDR < 0.05),]
        if(nrow(plotx) == 0){
          plotx <- data()
        }
        
        print("Plotx:")
        print(head(plotx))
        plotx <- plotx[which(plotx$FDR < 0.01),]
        plotx <- plotx[order(plotx$FDR, decreasing = F),]
        plotx <- plotx[order(plotx$Tau, decreasing = F),]
        if(nrow(plotx) > 20){
          plotx <- split(plotx, ifelse(plotx$Tau > 0, "Positive", "Negative"))
          plotx <- lapply(plotx, function(x){
            if(length(which(x$Tau > 0)) > 0){
              x <- x[order(x$Tau, decreasing = T),]
              x <- x[1:ifelse(nrow(x) > 10, 10, nrow(x)),]
            }else{
              x <- x[order(x$Tau, decreasing = F),]
              x <- x[1:ifelse(nrow(x) > 10, 10, nrow(x)),]
            }
          })
          plotx <- do.call(rbind.data.frame, plotx)
        }
        
        plotx$Drug_Name <- factor(plotx$Drug_Name, levels = unique(plotx[order(plotx$Tau, decreasing = F),"Drug_Name"]))
        plotx$FDR <- as.numeric(as.character(plotx$FDR))
        plotx$NegLog10FDR <- -log10(plotx$FDR)
        p <- NULL
        p <- ggplot(plotx, aes(NegLog10FDR, Tau, fill = Drug_Name, size = abs(Tau), label = Drug_Name))+
          geom_point(shape = 21, alpha = 0.8, color = "black")+
          facet_wrap(~Effect, ncol = 1, scales = "free_y")+
          scale_fill_manual(values = colorRampPalette(classcols)(length(unique(plotx$Drug_Name))), labels = ~stringr::str_wrap(.x, width = 1))+
          theme_classic(base_size = 12)+
          scale_size_continuous(range = c(1,10))+
          theme(strip.text.x = element_text(size = 14), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                legend.position = "none",
                panel.border = element_rect(colour = "black", fill = NA),
                axis.title = element_text(size = 15, face = "bold"))+
          ggtitle("Top 20 Drug Candidates\n(By Top 10 Postive/Negative Tau Scores)")+
          guides(size="none")
        return(p)
        
    })
    
    output$plot <- renderPlotly({
      PlotObject()
    })

    output$table <- renderDT({
      req(input$goi)
      data()[,which(colnames(data()) %in% c("Effect","LINCS_ID","Drug_Name","Mechanism","Drug_Status","Tau","ES","FDR","Drug_Specificity","Source_ID","PubChem_ID"))] %>% mutate_at(vars(Tau,FDR,ES,Drug_Specificity), funs(signif(., 3)))
    }, filter = "top", style="bootstrap", rownames = F, escape = FALSE, options = list(pageLength = 15))
    
    output$dltable <- downloadHandler(
      filename = function(){ paste("DRE_Drug_Repurposing_Result_",sys_time,".txt", sep = "")},
      content <- function(file){
        x <- data()[,which(colnames(data()) %in% c("Effect","LINCS_ID","Drug_Name","Mechanism","Drug_Status","Tau","ES","FDR","Drug_Specificity","Source_ID","PubChem_ID"))]
        write.table(x,file, row.names = F, quote = F, sep = "\t")
      })
    
    output$example <- renderImage({
      cfile <- list.files("DB/", pattern = "Online_Examples", full.names = T)
      list(src = cfile,
           contentType = 'image/png',
           width = "100%",
           height = "100%")
    }, deleteFile = FALSE)
    
}



