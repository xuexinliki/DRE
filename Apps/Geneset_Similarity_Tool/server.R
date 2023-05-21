function(input, output, session){
    corganism <- reactive({
      input$corganism
    })
    
    goi <- reactive({
      print("entering goi..")
      req(input$run)
      req(input$corganism)
      req(input$goi)
      print("goi..")
      x <- unlist(strsplit(input$goi, split = ","))
      x <- unlist(strsplit(x, split = "\n"))
      x <- unlist(strsplit(x, split = " "))
      x <- x[which(!is.na(x))]
      x <- x[which(x != "")]
      return(x)
    })
    
    sigsets <- reactive({
      print("entering sigsets..")
      showModal(modalDialog("Retrieving database..", footer=NULL))
    sigsets <- NULL
    sigsets <- readRDS(url(orgurl[which(orgurl$Species == gsub("_"," ",corganism())),"URL"], "rb"))
    sigsets <- lapply(sigsets, function(x){x <- toupper(x[!is.na(x)]); x <- unique(x)})
    removeModal()
    return(sigsets)
    })
    
    observeEvent(input$run, {
      print("entering observeEvent..")
      goi()
      sigsets()
      data()
    })
    
    data <- reactive({
      print("entering data..")
      req(corganism())

      showModal(modalDialog("Running analysis, the process may take a few minutes..", footer=NULL))
      
      cresult <- NULL
      genome.size <- length(unique(unlist(sigsets())))
      
      for(i in 1:length(sigsets())){
        print(i)
        coverlaps <- NULL
        coverlaps <- testGeneOverlap(newGeneOverlap(goi(), sigsets()[[i]], genome.size = genome.size))
        cresult <- rbind(cresult, data.frame(ID = names(sigsets())[i], P = coverlaps@pval, Similarity_Index = coverlaps@Jaccard))
      }
      
      cresult <- data.frame(MSigDB_Class = gsub(".*?\\|(.*?)\\|.*","\\1",cresult$ID),
                            Pathway = gsub(".*?\\|.*?\\|(.*)","\\1",cresult$ID),
                            cresult[,which(!colnames(cresult) %in% c("ID"))])
      
      cresult$Padj <- p.adjust(cresult$P, method = "BH")
      cresult <- cresult[order(cresult$Padj, decreasing = F),]
      print("Done!")
      removeModal()
      return(cresult)
    })
    
    
    PlotObject <- reactive({
      req(input$run)
      req(data())
      p <- NULL
        plotx <- data()
        plotx <- plotx[1:ifelse(nrow(plotx) > 20, 20, nrow(plotx)),]
        plotx$NegLogP <- -log(plotx$P)
        plotx$Pathway <- factor(plotx$Pathway, levels = rev(unique(plotx$Pathway)))
        p <- ggplot(plotx, aes(NegLogP, Pathway, fill = Similarity_Index, size = Similarity_Index))+
          geom_point(shape=21, stroke = 0.1)+theme_classic()+scale_fill_gradientn(colours = fillcols)
        return(p)
    })
    
    output$plot <- renderPlotly({
      PlotObject()
    })

    output$table <- renderDT({
      req(data())
      data() %>% mutate_at(vars(P,Padj,Similarity_Index), funs(signif(., 3)))
    }, filter = "top", style="bootstrap", rownames = F, escape = FALSE, options = list(pageLength = 15,
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



