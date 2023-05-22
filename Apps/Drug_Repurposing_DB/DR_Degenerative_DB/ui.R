fluidPage(titlePanel("DRE Drug-Repurposing for Disease Category: Degenerative Diseases"),
          sidebarLayout(
            sidebarPanel(
              selectizeInput('organism', 'Type to search for an organism*', choices = organism,selected = NULL,
                             options = list(placeholder = 'Please search and click to select an option',
                                            onInitialize = I('function() { this.setValue(""); }'))),
              selectizeInput("diseases","Type to search for a disease*", choices = diseases, selected = NULL, options = list(placeholder = 'Please search and click to select an option', onInitialize = I('function() { this.setValue(""); }'))),
              tags$i("*Press backspace to delete current search"),
              hr(),
h5("Pathway Category:"),
                             tags$a(href="https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2", "C2 (curated gene sets)", target="_blank"),
                             br(),
tags$a(href="https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C4", "C4 (computational gene sets)", target="_blank"),
br(),
                             tags$a(href="https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C5", "C5 (ontology gene sets)", target="_blank"),
                             br(),
                             tags$a(href="https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C6", "C6 (oncogenic signature gene sets)", target="_blank"),
                             br(),
                             tags$a(href="https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C7", "C7 (immunologic signature gene sets)", target="_blank"),
            # selectInput("pathway","Select a pathway", choices = pathway, selected = pathway[1])
h6("Note: The value of Tau score represents the extent of association of the drug with the pathway, such that a positive Tau (e.g. Tau = 100) of a pathway suggests that it is positively correlated with a drug and vice versa. Please note that both up/down regulated significant pathways (i.e., absolute Tau > 80) are shown for diseases, thus please refer to the Pathway_Description for further decripering of drug-pathway relationship."),
            ),
            mainPanel(
              tabsetPanel(
                tabPanel("Top Disease-Drug Associations",
                         br(),
                         # uiOutput("desc"),
                         useShinyjs(),
                         withSpinner(downloadablePlotlyUI('plot1', height = "1000px")),
                         # uiOutput("desc"),
                         # hr(),
                         # fluidRow(
                         #   conditionalPanel(condition = "input.type == 'Pathway'",
                         #                    withSpinner(downloadablePlotlyUI('plot2',height = "500px")))
                         # )
                         ),
                # downloadablePlotlyUI("plot", height = 600, width = "100%", inline = T)%>% withSpinner()),
                tabPanel("Query Table",br(),
                uiOutput("desc"),
                dataTableOutput("pathwaytable", width = "100%", height = 400)%>% withSpinner()),
                tabPanel("Example",br(),
                         div(withSpinner(imageOutput("example", width = "100%", height = "100%")), align = "center"))
              ))
          ),
)
