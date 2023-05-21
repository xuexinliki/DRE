function(input, output, session){
  example <- "DB/DRE_Drug_GSEA_Analysis_Example_Data_U2OS_PolyQ94_DEGs_Drug_Repurposing.txt"
  
  output$downloadData <- downloadHandler(
    filename = "DRE_Drug_GSEA_Analysis_Example_Data_U2OS_PolyQ94_DEGs_Drug_Repurposing.txt",
    content = function(file) {
      file.copy(example, file)
    },
    contentType = "text/txt"
  )
    
    crank <- reactive({
      req(input$run)
      req(input$cfile)
      cfile <- read.table(input$cfile$datapath[1], header = F, sep = "\t")
      cfile <- cfile[grep("Drug_Name|^Tau$", cfile[,1], ignore.case = T, invert = T),]
      crank <- NULL
      crank <- as.numeric(as.character(cfile[,2]))
      names(crank) <- toupper(cfile[,1])
      crank = sort(crank, decreasing = TRUE)
      return(crank)
    })
    
    observeEvent(input$run, {
      data()
    })
    
    data <- reactive({
      req(input$run)
      req(input$cfile)

      showModal(modalDialog("Running analysis which takes around a minute..", footer=NULL))
      cresult <- NULL
      cresult <- fgseaMultilevel(sigsets, stats = crank(), nPermSimple = 1000, eps = 0)
      cresult <- cresult[order(cresult$pval, decreasing = F),]
      cresult <- cresult[order(cresult$NES, decreasing = F),]
      removeModal()
      
      cpathways <- unique(cresult$pathway)
      updateSelectizeInput(session, inputId = 'cgsea', label = 'Select or type to search for a drug MOA family to view (press backspace to delete current selection)', choices = cpathways, selected = cpathways[1], server = T)
      
      return(cresult)
    })
    
    
    PlotObject <- reactive({
      req(input$cgsea)
      req(data())
      p <- NULL
        plotx <- data()
        plotx <- plotx[which(plotx$pathway == input$cgsea),]
        csets <- sigsets[input$cgsea]
        if(nrow(plotx) > 0){
          p <- plotEnrichment(sigsets[[input$cgsea]],crank(), gseaParam = 1)
          print(p)
        }
        return(p)
    })
    
    output$plot <- renderPlotly({
      PlotObject()
    })

    output$table <- renderDT({
      req(data())
      data() %>% mutate_at(vars(pval,padj,ES,NES), funs(signif(., 3)))
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



