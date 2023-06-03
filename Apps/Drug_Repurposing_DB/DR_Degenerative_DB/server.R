function(input, output, session){
  
  updateSelectizeInput(session, 'organism', 'Type to search for an organism*', choices = organism, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')))  
  data <- reactive({
    req(input$organism)
    
    print(paste("organism:", input$organism))
    plotx <- NULL
    db <- readRDS(paste("DB/GSEA_Result_Padj005_Tau_Mini/",gsub("\\s+","_",input$organism), "_GSEA_PAdj005_Tau_Abs80_Mini.RDS", sep = ""))
    db <- data.frame(Interaction = ifelse(db$Tau > 0, "Positive Tau", "Negative Tau"), db)
    db$NegLog10FDR <- -log10(db$FDR)
    db$Drug_Name <- paste(stringr::str_to_title(db$Drug_Name)," [Profile ID=",gsub("sig_","",db$Signature, ignore.case = T),"]", sep = "")
    db$Pathway <- gsub("'>","' target='_blank'>",db$Pathway, ignore.case = T)
    db$PubChem_ID <- gsub("'>","' target='_blank'>",db$PubChem_ID, ignore.case = T)
    db$LINCS_ID <- gsub("'>","' target='_blank'>",db$LINCS_ID, ignore.case = T)
    db$Pathway_PubMed <- gsub("'>","' target='_blank'>",db$Pathway_PubMed, ignore.case = T)
    db <- db[,which(!colnames(db) %in% c("Signature"))]
    db <- db[grepl("KEGG_", db$Pathway, ignore.case = T) | (db$Pathway_Category == "C2" & grepl("WP_", db$Pathway, ignore.case = T)),]
    db <- db[order(db$Tau, decreasing = F),]
    plotx$p1 <- db
    ccols <- colorRampPalette(classcols)(length(sort(unique(db$Drug_Name))))
    names(ccols) <- gsub(".*?>(.*?)<.*","\\1",sort(unique(db$Drug_Name)))
    plotx$ccols <- ccols
    cdiseases <- sort(unique(db$Disease_Type))
    plotx$cdiseases <- cdiseases
    cpathways <- sort(unique(db$Pathway_Category))
    plotx$cpathways <- cpathways
    # updateSelectizeInput(session, inputId = 'pathwayclass', choices = cpathways, server = TRUE)
    # updateSelectizeInput(session, inputId = 'diseasecats', choices = cdiseasecats, server = TRUE)
    updateSelectizeInput(session, "diseases","Type to search for a disease*", choices = cdiseases, selected = "", options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }')), server = TRUE)
    return(plotx)
  })
  
  
  pplot1 <- reactive({
    req(input$organism)
    req(data())
    req(input$diseases)
    print(paste("diseases:", input$diseases))
    x <- data()$p1[which(data()$p1$Disease_Type %in% input$diseases),]
    if(nrow(x)>0){
      x <- x[order(x$Tau, decreasing = F),]
      x <- x[which(x$FDR < 0.01),]
      x$Pathway <- gsub(".*>(.*)</a>$","\\1",x$Pathway, ignore.case = T)
      cgenes <- NULL
      cgenes <- split(x, ifelse(x$Tau > 0, "Positive", "Negative"))
      cgenes <- lapply(cgenes, function(x){
        if(length(which(x$Tau > 0)) > 0){
          x <- x[order(x$Tau, decreasing = T),]
          x <- x[1:ifelse(nrow(x) > 10, 10, nrow(x)),]
        }else{
          x <- x[order(x$Tau, decreasing = F),]
          x <- x[1:ifelse(nrow(x) > 10, 10, nrow(x)),]
        }
      })
      cgenes <- do.call(rbind.data.frame, cgenes)
      x <- x[which(x$Drug_Name %in% cgenes$Drug_Name),]
      x$Drug_Name <- factor(as.character(x$Drug_Name), levels = unique(cgenes$Drug_Name))
      cthreshold <- NULL
      cthreshold <- floor(min(abs(x$Tau)))
      x$Interaction <- paste(x$Interaction, " Pathways", sep = "")
      p <- NULL
      p <- ggplot(x, aes(Drug_Name, Tau, label = Pathway, fill = Drug_Name))+
        geom_point(alpha = 0.9, aes(size = NegLog10FDR))+xlab("")+ylab("Tau Score")+
        facet_wrap(~Interaction, ncol = 1, scales = "free_y")+
        scale_fill_manual(values = data()$ccols, labels = ~stringr::str_wrap(.x, width = 1))+
        scale_shape_manual(values = c(Approved = 21,Experimental=23,Withdrawn=24), labels = ~ stringr::str_wrap(.x, width = 30))+
        scale_size_manual(labels = ~ stringr::str_wrap(.x, width = 30))+
        theme_classic(base_size = 12)+
        scale_size_continuous(range = c(1,10))+
        theme(strip.text.x = element_text(size = 14), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              legend.position = "none",
              panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
        # geom_hline(yintercept=0, linetype="dashed", color = "black")+
        guides(color=guide_legend(ncol=2), size="none", shape="none")
      # guide_legend(title="-log10(FDR)"), guide_legend(title="Drug Status")
      ggplotly(p)  %>%
        config(
          toImageButtonOptions = list(
            format = c("png")))
    }
  })
  
  corganism <- reactive({
    req(input$organism)
    input$organism
  })
  
  callModule(downloadablePlotly, 
             id = 'plot1', 
             plot = pplot1, 
             filename = function(name) {paste0("OFCDE_Barplot_Input_",input$diseases,"_",corganism(),".csv")}, 
             content = function(file) {write.csv(data()$p1, file, quote = F)})
  
  output$pathwaytable <- renderDataTable(
    data()$p1[which(data()$p1$Disease_Type %in% input$diseases),] %>% mutate_at(vars(Tau,ES,FDR,Drug_Specificity), funs(signif(., 3))),
    filter = "top",
    extensions = 'Buttons',
    style="bootstrap",
    class = "compact",
    rownames = F,
    escape = FALSE,
    options = list(pageLength = 10,
                   buttons = c('copy', 'excel'),
                   dom = 'Blfrtip')) # <lf<\"datatables-scroll\"t>ipr>
  
  output$desc <- renderUI(HTML(markdown::renderMarkdown(text = "Note: The value of Tau score represents the extent of association of the drug with the pathway, such that a positive Tau (e.g. Tau = 100) of a pathway suggests that it is positively correlated with a drug and vice versa. Please note that both up/down regulated significant pathways (i.e., absolute Tau > 80) are shown for diseases, thus please refer to the Pathway_Description for further decripering of drug-pathway relationship.")))
  
  output$example <- renderImage({
    cfile <- list.files("DB/", pattern = "Online_Examples", full.names = T)
    list(src = cfile,
         contentType = 'image/png',
         width = "100%",
         height = "100%")
  }, deleteFile = FALSE)
  
}