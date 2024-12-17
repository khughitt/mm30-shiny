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
        "OS",
        fluidRow(
          selectInput("surv_os_feat_type", "Genes / Pathways:", 
                      c("Genes"="genes", "Pathways"="pathways"))
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
                  width=3, 
                  uiOutput("surv_os_gene_summary_html"),
                )
              ),
              withSpinner(plotOutput("surv_os_plot", height="800px")),
              tagList(
                tags$hr(),
                tags$h3("Survival regression results:")
              ),
              tableOutput("surv_os_gene_details_tbl"),
              tagList(
                tags$hr(),
                tags$h3("TCGA survival associations:")
              ),
              tableOutput("surv_os_gene_tcga_tbl")
            ),
          ),
        ),
        #---------------------------------------
        # Overall survival > pathways
        #---------------------------------------
        conditionalPanel(
          condition = "input.surv_os_feat_type == 'pathways'",
          fluidRow(
            column(
              width=6,
              withSpinner(DTOutput("surv_os_pathway_tbl"))
            ),
            column(
              width=6,
              fluidRow(
                column(
                  width=3, 
                  uiOutput("surv_os_pathway_summary_html"),
                )
              ),
              #withSpinner(plotOutput("surv_os_plot", height="800px")),
              tagList(
                tags$hr(),
                tags$h3("Survival regression results:")
              ),
              tableOutput("surv_os_pathway_details_tbl"),
              tagList(
                tags$hr(),
                tags$h3("Genes:"),
              ),
              tableOutput("surv_os_pathway_gene_tbl"),
            ),
          ),
        ),
      ),
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
      tabPanel(
        "Settings",
        fluidRow(
        )
      )
    )
  )
}
