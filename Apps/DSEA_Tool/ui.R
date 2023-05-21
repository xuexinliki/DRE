fluidPage(titlePanel("DRE: Drug-Set Enrichment Analysis"),
          sidebarLayout(
            sidebarPanel(
              fileInput('cfile', 'Upload a two column (no header) tab-delimited .txt file containing drug names with profile ID (column 1) and their corresponding Tau scores(column 2) from a post-drug repurposing analysis output.',
                        accept=c('text/txt', 'text/txt')),
              downloadButton("downloadData", label = "Download an Example Input Data"),
              hr(),
              h5("Example of U2OS PolyQ94 DEGs Drug Repurposing Analysis output:"),
              h6("Ly-341495 [Profile ID=74]	-99.811148339933"),
              h6("Rivaroxaban [Profile ID=6584]	-99.8081023454158"),
              h6(".."),
              h6("Rifampicin [Profile ID=20267]	-99.7410904660372"),
              actionButton("run", "Start Analysis"),
            ),
            mainPanel(
              fluidRow(column(12,selectInput(inputId = 'cgsea', label = 'Select or type to search for a pathway to view (press backspace to delete current selection)',width="100%",
                                             choices = "", multiple = F, selected = "", selectize = T))),
              tabsetPanel(
                tabPanel("Plot",
                         br(),
                         fluidRow(column(12,withSpinner(plotlyOutput('plot', width = "100%", height = 500))))
                ),
                tabPanel("Table",br(),
                         # hr(),
                         fluidRow(column(12, DTOutput('table')))),
                tabPanel("Example",br(),
                         div(withSpinner(imageOutput("example", width = "100%", height = "100%")), align = "center"))
                
              ))),
)
