function(input, output, session){
  
  updateSelectizeInput(session, "organism","Type to search for an organism*", choices = organism, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')))
  updateSelectizeInput(session, "type","Select a data type", choices = c("Pathway","Drug"), options = list(placeholder = 'Please click to select an option', onInitialize = I('function() { this.setValue(""); }')))

  data <- reactiveValues(db = NULL, cmap = NULL, sigsets = NULL, gsearesult = NULL)
  
  observeEvent(input$organism, {
    if(!is.null(input$organism)){
      if(input$organism != ""){
        showModal(modalDialog("Reading in database..", footer=NULL))
        x <- NULL
        x <- readRDS(url(orgurl[which(orgurl$Organism == input$organism),"URL"],"rb"))
        x$Drug_Name <- paste(stringr::str_to_title(x$Drug_Name)," [Profile ID=",gsub("sig_","",x$Signature, ignore.case = T),"]", sep = "")
        x <- x[,which(!colnames(x) %in% c("Signature"))]
        colnames(x)[which(colnames(x) == "MSigDB_Class")] <- "Pathway_Category"
        x[which(!is.na(x$PubChem_ID)),"PubChem_ID"] <- paste0("<a href='https://pubchem.ncbi.nlm.nih.gov/compound/",x[which(!is.na(x$PubChem_ID)),"PubChem_ID"],"' target='_blank'>",x[which(!is.na(x$PubChem_ID)),"PubChem_ID"],"</a>")
        x[which(!is.na(x$LINCS_ID)),"LINCS_ID"] <- paste0("<a href='http://www.ilincs.org/ilincs/perturbagen/compound/",x[which(!is.na(x$LINCS_ID)),"LINCS_ID"],"' target='_blank'>",x[which(!is.na(x$LINCS_ID)),"LINCS_ID"],"</a>")
        cref <- csummary[which(csummary$Species == gsub("_"," ",input$organism)),]
        x$Pathway_Description <- unlist(cref[match(x$Pathway, cref$GeneSet_Name),"GeneSet_Description"])
        x$Pathway_PubMed <- paste("<a href='",unlist(cref[match(x$Pathway, cref$GeneSet_Name),"PubMedURL"]),"' target='_blank'>",unlist(cref[match(x$Pathway, cref$GeneSet_Name),"PubMed_ID"]),"</a>", sep = "")
        x$Pathway_PubMed[which(x$Pathway_PubMed == "<a href='https://pubmed.ncbi.nlm.nih.gov/' target='_blank'></a>")] <- NA
        x$Pathway_Name <- x$Pathway
        x$Pathway <- paste("<a href='",unlist(cref[match(x$Pathway, cref$GeneSet_Name),"MSigdbURL"]),"' target='_blank'>",x$Pathway,"</a>", sep = "")
        x <- data.frame(Effect = ifelse(x$Tau > 0, "Postive Tau", "Negative Tau"), x)
        x$Effect[which(x$Tau == 0)] <- "Zero Tau"
        x <- x[order(x$Tau, decreasing = F),]
        data$db <- x
        sigsets <- NULL
        clink <- NULL
        clink <- reflinks[match(gsub(" ","_",input$organism, ignore.case = T), gsub(" ","_",reflinks$Species, ignore.case = T)),"URL"]
        clink <- gsub("dl=0","dl=1",clink)
        sigsets <- readRDS(url(clink,"rb"))
        sigsets <- lapply(sigsets, function(x){x <- toupper(x[!is.na(x)]); x <- unique(x)})
        data$sigsets <- sigsets
        removeModal()
      }
    }
  })
  
observeEvent(input$type, {
    if(!is.null(data$db)){
      if(input$type != ""){
        plotx <- NULL
        if(input$type == "Drug"){
          cdrugnames <- sort(unique(data$db$Drug_Name))
          updateSelectizeInput(session, "drugname","Type to search for a drug name*", choices =  cdrugnames, selected = input$drugname, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')), server = TRUE)
          showTab(inputId = "tabs", target = "GSEA Plot")
        }else if(input$type == "Pathway"){
          cpathwaycats <- sort(unique(data$db$Pathway_Category))
          updateSelectizeInput(session, "pathwayclass2","Type to search for a pathway category*", choices = cpathwaycats, selected = input$pathwayclass2, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')), server = TRUE)
          hideTab(inputId = "tabs", target = "GSEA Plot")
        }
      }
    }
  })
  
observeEvent(input$drugname, {
  if(!is.null(data$db)){
    if(input$drugname != ""){
      plotx <- NULL
      if(input$type == "Drug"){
        cpathwaycats <- sort(unique(data$db[which(data$db$Drug_Name == input$drugname),"Pathway_Category"]))
        updateSelectizeInput(session, "pathwayclass1","Type to search for a pathway category*", choices = cpathwaycats, selected = input$pathwayclass1, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')), server = TRUE)
      }
    }
  }
})

observeEvent(input$pathwayclass1, {
  if(!is.null(data$db)){
    print("input$type,input$pathwayclass1,input$drugname:")
    print(paste(input$type,input$pathwayclass1,input$drugname, sep = ","))
    if(input$type == "Drug" & input$pathwayclass1 != "" & input$drugname != ""){
        cpathways <- sort(unique(data$db[which(data$db$Drug_Name == input$drugname & data$db$Pathway_Category == input$pathwayclass1),"Pathway_Name"]))
        if(length(cpathways) > 0){
          updateSelectizeInput(session, "pathway","Type to search for a pathway to view its GSEA plot*", choices = cpathways, selected = input$pathway, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')), server = TRUE)
        }
    }
  }
})

observeEvent(input$pathwayclass2, {
  if(!is.null(data$db)){
    print("input$type,input$pathwayclass2,input$pathway2:")
    print(paste(input$type,input$pathwayclass2, sep = ","))
    if(input$type == "Pathway" & input$pathwayclass2 != ""){
      cpathways <- sort(unique(data$db[which(data$db$Pathway_Category == input$pathwayclass2),"Pathway_Name"]))
      if(length(cpathways) > 0){
        updateSelectizeInput(session, "pathway2","Type to search for a pathway*", choices = cpathways, selected = input$pathway2, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')), server = TRUE)
      }
    }
  }
})

  pplot1 <- reactive({
    req(input$organism)
    req(input$type)
    req(!is.null(data$db))
    req((input$type == "Drug" & input$drugname != "" & input$pathwayclass1 != "") | (input$type == "Pathway" & !input$pathwayclass2 %in% c("",NULL) & input$pathway2 != ""))
    print(paste("input$pathwayclass2:", input$pathwayclass2))
    print(paste("type:", input$type))
    if(input$type == "Drug"){
      print(paste("drugname:", input$drugname))
       if(input$drugname != ""){
         x <- NULL
         x <- data$db[which(data$db$Drug_Name == input$drugname & data$db$Pathway_Category == input$pathwayclass1),]
         print("nrow(x):")
         print(nrow(x))
         plotx <- NULL
         if(nrow(x) > 0){
           x$NegLog10FDR <- -log10(x$FDR)
           x$Effect <- factor(x$Effect, levels = c("Postive Tau", "Negative Tau", "Zero Tau"))
           # ccols <- NULL
           # ccols <- colorRampPalette(classcols)(length(sort(unique(x$Mechanism))))
           # names(ccols) <- gsub(".*?>(.*?)<.*","\\1",sort(unique(x$Mechanism)))
           p <- NULL
           p <- ggplot(x, aes(NegLog10FDR, Tau, fill = Pathway, size = abs(Tau), label = Pathway))+
             geom_point(shape = 21, alpha = 0.8, color = "black")+
             facet_wrap(~Effect, ncol = 1, scales = "free_y")+
             scale_fill_manual(values = colorRampPalette(classcols)(length(unique(x$Pathway))), labels = ~stringr::str_wrap(.x, width = 1))+
             theme_classic(base_size = 12)+
             scale_size_continuous(range = c(1,10))+
             theme(strip.text.x = element_text(size = 14), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   strip.background = element_blank(),
                   legend.position = "none",
                   panel.border = element_rect(colour = "black", fill = NA),
                   axis.title = element_text(size = 15, face = "bold"))+
             ggtitle("Drug-Associated Pathway Signatures")+
             guides(size="none")
          
           ggplotly(p)  %>% 
             config(
               toImageButtonOptions = list(
                 format = c("png")))
       }

      }else{
        return(NULL)
        # ggplot() +                      # Draw ggplot2 plot with text only
        #   annotate("text",
        #            x = 3,
        #            y = 3,
        #            size = 4,
        #            label = "There is no associated pathway for this drug at the current\nsignificance level. Please select another drug.") + 
        #   theme_minimal()
      }
    }else if(input$type == "Pathway"){
      print(paste("pathwayclass:", input$pathwayclass2))
      if(input$pathwayclass2 != "" & input$pathway2 != ""){
      print(paste("current pathwayclass:", input$pathwayclass2))
      print(paste("current pathway:", input$pathway2))
      x <- data$db[which(data$db$Pathway_Category %in% input$pathwayclass2 & data$db$Pathway_Name == input$pathway2),]
      x <- x[order(x$Tau, decreasing = F),]
      x <- split(x, ifelse(x$Tau > 0, "Positive", "Negative"))
      x <- lapply(x, function(x){
        if(length(which(x$Tau > 0)) > 0){
          x <- x[order(x$Tau, decreasing = T),]
          x <- x[1:ifelse(nrow(x) > 10, 10, nrow(x)),]
        }else{
          x <- x[order(x$Tau, decreasing = F),]
          x <- x[1:ifelse(nrow(x) > 10, 10, nrow(x)),]
        }
      })
      x <- do.call(rbind.data.frame, x)
      
      if(nrow(x) > 0){
        x$NegLog10FDR <- -log10(x$FDR)
        x$Effect <- factor(x$Effect, levels = c("Postive Tau", "Negative Tau", "Zero Tau"))
        x <- x[order(x$NegLog10FDR, decreasing = F),]
        x$Drug_Name <- factor(x$Drug_Name, levels = unique(x$Drug_Name))
        
        # ccols <- NULL
        # ccols <- colorRampPalette(classcols)(length(sort(unique(x$Mechanism))))
        # names(ccols) <- gsub(".*?>(.*?)<.*","\\1",sort(unique(x$Mechanism)))
        p <- NULL
        p <- ggplot(x, aes(NegLog10FDR, Tau, fill = Drug_Name, size = abs(Tau), label = Drug_Name))+
          geom_point(shape = 21, alpha = 0.8, color = "black")+
          facet_wrap(~Effect, ncol = 1, scales = "free_y")+
          scale_fill_manual(values = colorRampPalette(classcols)(length(unique(x$Drug_Name))), labels = ~stringr::str_wrap(.x, width = 1))+
          theme_classic(base_size = 12)+
          scale_size_continuous(range = c(1,10))+
          theme(strip.text.x = element_text(size = 14), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                legend.position = "none",
                panel.border = element_rect(colour = "black", fill = NA),
                axis.title = element_text(size = 15, face = "bold"))+
          ggtitle("Top 20 Drugs with Positive/Negative Tau")+
          guides(size="none")
        
        ggplotly(p)  %>% 
          config(
            toImageButtonOptions = list(
              format = c("png")))
      }else{
        return(NULL)
        # ggplot() +                      # Draw ggplot2 plot with text only
        #   annotate("text",
        #            x = 3,
        #            y = 3,
        #            size = 4,
        #            label = "There is no associated pathway for this drug at the current\nsignificance level. Please select another drug.") + 
        #   theme_minimal()
      }
      # cpaths <- NULL
      # cpaths <- unique(x$Pathway)
      # cpaths <- cpaths[1:ifelse(length(cpaths) > 20, 20, length(cpaths))]
      # x <- x[which(x$Drug_Name %in% cdrugs$Drug_Name & x$Pathway %in% cdrugs$Pathway),]
      # x <- x[,c("Drug_Name","Pathway","Mechanism","Drug_Status","Group")]
      # x <- x %>%
      #   make_long(Drug_Status,Mechanism,Drug_Name,Effect,Pathway)
      # colnames(x) <- c("Level_I","Stage_I","Level_II","Stage_II")
      # 
      # ccols <- NULL
      # ccols <- colorRampPalette(groupcols)(100)
      # x$Stage_I <- factor(x$Stage_I)
      # p <- ggplot(x, aes(x = Level_I, 
      #                    next_x = Level_II, 
      #                    node = Stage_I, 
      #                    next_node = Stage_II,
      #                    fill = Stage_I)) + # label = Level_I
      #   geom_sankey(flow.alpha = 0.8, node.color = 1) +
      #   theme_sankey(base_size = 20)+
      #   theme(legend.position = "none")+ # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      #   # geom_sankey_label(size = 3.5, color = 1, fill = "white") +
      #   scale_fill_manual(values = ccols) +
      #   theme(legend.position = "none")+
      #   xlab("")+ggtitle("Top 20 Drug-Pathway Associations")
      
      # p %>%
      #   ggplotly() %>%
      #   config(
      #     toImageButtonOptions = list(
      #       format = c("png")))
      }
    }
  })
  
  callModule(downloadablePlotly, 
             id = 'plot1', 
             plot = pplot1, 
             filename = function(name) {paste0("DRE_Barplot_Input_",input$drugname,"_",input$organism,".csv")}, 
             content = function(file) {write.csv(data$db, file, quote = F)})
  
  pplot2 <- reactive({
    req(input$organism)
    req(input$type)
    req(!is.null(data$db))
    req(input$type == "Drug" & input$drugname != "" & input$pathwayclass1 != "" & input$pathway != "")
      print(paste("pathway:", input$pathway))
      if(input$pathway != ""){
        print("setting ccol..")
        ccol <- which(colnames(cmap) %in% paste("sig_",ref$sig_id, sep = "")[which(ref$Drug_Name == input$drugname)])
        crank <- NULL
        crank <- cmap[,ccol]
        names(crank) <- row.names(cmap)
        crank = sort(crank, decreasing = TRUE)
        csets <- data$sigsets[which(gsub(".*?\\|(.*?)\\|.*","\\1",names(data$sigsets)) %in% input$pathwayclass1)]
        names(csets) <- gsub(".*?\\|.*?\\|(.*)","\\1",names(csets))
        csets <- csets[which(names(csets) %in% input$pathway)]
        cresult <- fgsea::fgseaMultilevel(csets, stats = crank, nPermSimple = 10000, eps = 0)
        cresult <- data.frame(cresult)
        cresult <- cresult[,which(!colnames(cresult) %in% c("leadingEdge","size"))]
        cresult <- cresult[order(cresult$padj, decreasing = F),]
        data$gsearesult <- cresult
        print("cresult:")
        print(cresult)
        p <- plotEnrichment(csets[[input$pathway]],crank)
        return(p)
      }else{
        return(NULL)
        # ggplot() +                      # Draw ggplot2 plot with text only
        #   annotate("text",
        #            x = 3,
        #            y = 3,
        #            size = 4,
        #            label = "There is no associated pathway for this drug at the current\nsignificance level. Please select another drug.") + 
        #   theme_minimal()
    }
  })
  
  callModule(downloadablePlotly, 
             id = 'plot2', 
             plot = pplot2, 
             filename = function(name) {paste0("DRE_GSEA_Plot_Input_",input$drugname,"_",input$organism,".csv")}, 
             content = function(file) {write.csv(data$db, file, quote = F)})
  
  output$drugtable <- renderDataTable(
  data$db[which(data$db$Drug_Name == input$drugname & data$db$Pathway_Category == input$pathwayclass1),which(colnames(data$db) %in% c("Pathway_Category","Pathway","Pathway_Description","Tau","ES","FDR","Pathway_PubMed","Compound_ID","LINCS_ID"))] %>% mutate_at(vars(Tau,ES,FDR), funs(signif(., 3))),
                                  filter = "top",
                                  extensions = 'Buttons',
                                  style="bootstrap",
                                  class = "compact",
                                  rownames = F,
                                  escape = FALSE,
                                  options = list(pageLength = 10,
                                                 buttons = c('copy', 'excel'),
                                                 dom = 'Blfrtip')) # <lf<\"datatables-scroll\"t>ipr>
  
  
  output$pathwaytable <- renderDataTable(
    data$db[which(data$db$Pathway_Category == input$pathwayclass2),which(colnames(data$db) %in% c("Drug_Name","Mechanism","Drug_Status","Pathway_Category","Pathway","Pathway_Description","Tau","ES","FDR","Drug_Specificity","Pathway_PubMed","Compound_ID","LINCS_ID"))] %>% mutate_at(vars(Tau,ES,FDR,Drug_Specificity), funs(signif(., 3))),
    filter = "top",
    extensions = 'Buttons',
    style="bootstrap",
    class = "compact",
    rownames = F,
    escape = FALSE,
    options = list(pageLength = 10,
                   buttons = c('copy', 'excel'),
                   dom = 'Blfrtip')) # <lf<\"datatables-scroll\"t>ipr>
  
  output$drugdesc <- renderUI(HTML(markdown::renderMarkdown(text = paste("- Drug Name:",input$drugname,"\n- Mechanism:",unique(data$db[which(data$db$Drug_Name == input$drugname),"Mechanism"]),"\n- Drug Status:", unique(data$db[which(data$db$Drug_Name == input$drugname),"Drug_Status"]), "\n - Drug Specificity:",unique(data$db[which(data$db$Drug_Name == input$drugname),"Drug_Specificity"]),"\n- PubChem ID:", unique(data$db[which(data$db$Drug_Name == input$drugname),"PubChem_ID"])))))

  # output$pathwaydesc <- renderUI(HTML(markdown::renderMarkdown(text = paste("- Pathway:",input$pathway,"- MSigDB Class:",unique(data$db[which(data$db$Pathway == input$pathway),"Pathway_Category"])))))
  
  output$gseatable <- renderDT({
    req(data$gsearesult)
    print("nrow(data$gsearesult):")
    print(nrow(data$gsearesult))
    data$gsearesult %>% mutate_at(vars(pval,padj,ES,NES,log2err), funs(signif(., 3)))
  }, filter = "top", style="bootstrap", rownames = F, escape = FALSE, options = list(pageLength = 10,
                                                                                     deferRender = TRUE,
                                                                                     scrollY = 400,
                                                                                     scrollX = TRUE,
                                                                                     scroller = TRUE,
                                                                                     autoWidth = TRUE))
  
  output$example <- renderImage({
    cfile <- list.files("DB/", pattern = "Online_Examples", full.names = T)
    list(src = cfile,
         contentType = 'image/png',
         width = "100%",
         height = "100%")
  }, deleteFile = FALSE)
  
}



