ui <- function(request) {
  tagList(
    tags$head(
      includeCSS("resources/styles.css"),
      includeCSS("https://fonts.googleapis.com/css?family=Roboto+Mono&display=swap"),
      tags$style(
        HTML(".navbar-brand { font-family: 'Roboto Mono', monospace; }")
      )
    ),

    navbarPage(
      id="tab",
      theme=shinytheme("flatly"),
      title=textOutput("page_title"),
      windowTitle="MM30",
      ######################################## 
      # Overall survival
      ######################################## 
      tabPanel(
        "Survival (OS)",
        fluidRow(
          selectInput("surv_os_feat_type", "Genes / Gene sets:", 
                      c("Genes"="genes", "Gene sets"="gene_sets"))
        ),
        #---------------------------------------
        # Overall survival > genes
        #---------------------------------------
        conditionalPanel(
          condition = "input.surv_os_feat_type == 'genes'",
          fluidRow(
            column(
              width=6,
              withSpinner(DTOutput("surv_os_gene_scores_tbl"))
            ),
            column(
              width=6,
              fluidRow(
                column(
                  width=6, 
                  uiOutput("surv_os_gene_summary_html"),
                )
              ),
              withSpinner(plotOutput("surv_os_gene_plot", height="800px")),
              withSpinner(
                tagList(
                  tags$hr(),
                  tags$h3("Survival regression results:"),
                  tableOutput("surv_os_gene_details_tbl"),
                  tags$hr(),
                  tags$h3("TCGA survival associations:"),
                  tableOutput("surv_os_gene_tcga_tbl")
                ),
              ),
            ),
          ),
        ),
        #---------------------------------------
        # Overall survival > gene sets
        #---------------------------------------
        conditionalPanel(
          condition = "input.surv_os_feat_type == 'gene_sets'",
          fluidRow(
            column(
              width=6,
              withSpinner(DTOutput("surv_os_gene_set_scores_tbl"))
            ),
            column(
              width=6,
              fluidRow(
                column(
                  width=6, 
                  uiOutput("surv_os_gene_set_summary_html"),
                )
              ),
              withSpinner(plotOutput("surv_os_gene_set_plot", height="800px")),
              withSpinner(
                tagList(
                  tags$hr(),
                  tags$h3("Survival regression results:"),
                  tableOutput("surv_os_gene_set_details_tbl"),
                  tags$hr(),
                  tags$h3("Gene set members:"),
                  tableOutput("surv_os_gene_set_gene_tbl"),
                ),
              )
            ),
          ),
        ),
      ),
      ######################################## 
      # Disease stage
      ######################################## 
      tabPanel(
        "Disease Stage",
        fluidRow(
          selectInput("stage_feat_type", "Genes / Gene sets:", 
                      c("Genes"="genes", "Gene sets"="gene_sets"))
        ),
        #---------------------------------------
        # Disease stage > genes
        #---------------------------------------
        conditionalPanel(
          condition = "input.stage_feat_type == 'genes'",
          fluidRow(
            column(
              width=6,
              withSpinner(DTOutput("stage_gene_scores_tbl"))
            ),
            column(
              width=6,
              fluidRow(
                column(
                  width=3, 
                  uiOutput("stage_gene_summary_html"),
                )
              ),
              withSpinner(plotOutput("stage_gene_plot", height="800px")),
              withSpinner(
                tagList(
                  tags$hr(),
                  tags$h3("Logistic regression results:"),
                  tableOutput("stage_gene_details_tbl"),
                  tags$hr(),
                  tags$h3("TCGA survival associations:"),
                  tableOutput("stage_gene_tcga_tbl")
                ),
              ),
            ),
          ),
        ),
        #---------------------------------------
        # Disease stage > gene sets
        #---------------------------------------
        conditionalPanel(
          condition = "input.stage_feat_type == 'gene_sets'",
          fluidRow(
            column(
              width=6,
              withSpinner(DTOutput("stage_gene_set_scores_tbl"))
            ),
            column(
              width=6,
              fluidRow(
                column(
                  width=6, 
                  uiOutput("stage_gene_set_summary_html"),
                )
              ),
              withSpinner(plotOutput("stage_gene_set_plot", height="800px")),
              withSpinner(
                tagList(
                  tags$hr(),
                  tags$h3("Regression results:"),
                  tableOutput("stage_gene_set_details_tbl"),
                  tags$hr(),
                  tags$h3("Genes:"),
                  tableOutput("stage_gene_set_gene_tbl"),
                ),
              )
            ),
          ),
        ),
      ),
      ######################################## 
      # Co-expression
      ######################################## 
      tabPanel(
        "Co-expression",
        fluidRow(
          selectInput("coex_feat_type", "Genes / Gene sets:", 
                      c("Genes"="genes", "Gene sets"="gene_sets"))
        ),
        #---------------------------------------
        # co-expression > genes
        #---------------------------------------
        conditionalPanel(
          condition = "input.coex_feat_type == 'genes'",
          fluidRow(
            column(
              width=3,
              fluidRow(
                selectizeInput("coex_gene1", "Gene A:", choices=gene_coex_opts, selected=gene_coex_opts[1]),
                selectizeInput("coex_gene2", "Gene B:", choices=gene_coex_opts, selected=gene_coex_opts[2])
              ),
            ),
            column(
              width=9,
              withSpinner(plotOutput("gene_coex_plot", width="1200px", height="1000px"))
            )
          )
        )
      ),
      ######################################## 
      # Cell Lines
      ######################################## 
      tabPanel(
        "Cell Lines",
        fluidRow(
          column(
              width=6,
              withSpinner(DTOutput("hmcl_expr_tbl"))
          ),
          column(
            width=6,
            withSpinner(plotOutput("hmcl_hist_plot", height="800px")),
            fluidRow(
              column(
                width=3, 
                #uiOutput("hmcl_gene_summary_html"),
              ),
              column(
                width=3,
              ),
            ),
          ),
        ),
      ),
      ######################################## 
      # Settings
      ######################################## 
      # tabPanel(
      #   "Settings",
      #   fluidRow(
      #   )
      # )
    )
  )
}
