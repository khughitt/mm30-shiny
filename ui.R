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
              uiOutput("surv_os_gene_summary_html"),
              withSpinner(plotOutput("surv_os_gene_plot", height="800px")),
              tagList(
                tags$hr(),
                column(
                  width=6, 
                  tags$h3("Survival regression results:"),
                  tableOutput("surv_os_gene_details_tbl"),
                ),
                column(
                  width=6, 
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
              uiOutput("surv_os_gene_set_summary_html"),
              withSpinner(plotOutput("surv_os_gene_set_plot", height="800px")),
              tagList(
                tags$hr(),
                tags$h3("Survival regression results:"),
                tableOutput("surv_os_gene_set_details_tbl"),
                tags$hr(),
                tags$h3("Gene set members:"),
                tableOutput("surv_os_gene_set_gene_tbl"),
              ),
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
              withSpinner(DTOutput("gene_stage_scores_tbl"))
            ),
            column(
              width=6,
              uiOutput("gene_stage_summary_html"),
              withSpinner(plotOutput("gene_stage_plot", height="800px")),
              tagList(
                tags$hr(),
                column(
                  width=6, 
                  tags$h3("Logistic regression results:"),
                  tableOutput("gene_stage_details_tbl"),
                ),
                column(
                  width=6, 
                  tags$h3("TCGA survival associations:"),
                  tableOutput("gene_stage_tcga_tbl")
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
              withSpinner(DTOutput("gene_set_stage_scores_tbl"))
            ),
            column(
              width=6,
              uiOutput("gene_set_stage_summary_html"),
              withSpinner(plotOutput("gene_set_stage_plot", height="800px")),
              tagList(
                tags$hr(),
                tags$h3("Regression results:"),
                tableOutput("gene_set_stage_details_tbl"),
                tags$hr(),
                tags$h3("Genes:"),
                tableOutput("gene_stage_set_gene_tbl"),
              ),
            ),
          ),
        ),
      ),
      ######################################## 
      # Disease stage transitions
      ######################################## 
      tabPanel(
        "Stage Transitions",
        fluidRow(
            column(
              width=2,
              selectInput("transition_feat_type", "Genes / Gene sets:", 
                          c("Genes"="genes", "Gene sets"="gene_sets"))
            ),
            column(
              width=2,
              selectInput("transition_opt", "Transition:", 
                          c("Healthy vs. MGUS"="Healthy_MGUS",
                            "MGUS vs. SMM"="MGUS_SMM",
                            "SMM vs. MM"="SMM_MM",
                            "MM vs. RRMM"="MM_RRMM",
                            "Healthy/MGUS/SMM vs. MM/RRMM"="early_vs_late", 
                            "Healthy/MGUS vs. SMM/MM/RRMM"="before_vs_after_smm",
                            "Pre- vs. Post-relapse"="before_vs_after_relapse"))
            )
        ),
        #---------------------------------------
        # Disease stage transitions > genes
        #---------------------------------------
        conditionalPanel(
          condition = "input.transition_feat_type == 'genes'",
          fluidRow(
            column(
              width=4,
              withSpinner(DTOutput("gene_stage_transitions_tbl"))
            ),
            column(
              width=6,
              withSpinner(plotOutput("gene_stage_transitions_plot", height="800px")),
            ),
            column(
              width=2,
              tags$h3("Individual datasets:"),
              tableOutput("gene_stage_transitions_details_tbl")
            )
          )
        )
      ),
      ######################################## 
      # Treatment
      ######################################## 
      tabPanel(
        "Treatment",
        fluidRow(
          selectInput("treatment_feat_type", "Genes / Gene sets:", 
                      c("Genes"="genes", "Gene sets"="gene_sets"))
        ),
        #---------------------------------------
        # Treatment > genes
        #---------------------------------------
        conditionalPanel(
          condition = "input.treatment_feat_type == 'genes'",
          fluidRow(
            column(
              width=6,
              withSpinner(DTOutput("treatment_gene_scores_tbl"))
            ),
            column(
              width=6,
              uiOutput("treatment_gene_summary_html"),
              withSpinner(plotlyOutput("treatment_gene_plot", height="800px")),
              tagList(
                tags$hr(),
                column(
                  width=6, 
                  tags$h3("Treatment association test results:"),
                  tableOutput("treatment_gene_details_tbl"),
                ),
                column(
                  width=6, 
                  tags$h3("TCGA survival associations:"),
                  tableOutput("treatment_gene_tcga_tbl")
                ),
              ),
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
              width=2,
              fluidRow(
                selectizeInput("coex_gene1", "Gene A:", choices=gene_coex_opts, selected=gene_coex_opts[1]),
                selectizeInput("coex_gene2", "Gene B:", choices=gene_coex_opts, selected=gene_coex_opts[2])
              ),
            ),
            column(
              width=9,
              withSpinner(plotOutput("gene_coex_plot", width="2400px", height="1000px"))
            )
          )
        ),
        #---------------------------------------
        # co-expression > gene sets
        #---------------------------------------
        conditionalPanel(
          condition = "input.coex_feat_type == 'gene_sets'",
          fluidRow(
            column(
              width=3,
              fluidRow(
                selectizeInput("coex_gene_set1", "Gene Set A:", choices=gene_set_coex_opts, selected=gene_set_coex_opts[1]),
                selectizeInput("coex_gene_set2", "Gene Set B:", choices=gene_set_coex_opts, selected=gene_set_coex_opts[2])
              ),
            ),
            column(
              width=9,
              withSpinner(plotOutput("gene_set_coex_plot", width="1200px", height="1000px"))
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
      # Datasets
      ######################################## 
      tabPanel(
        "Datasets",
        fluidRow(
          column(
              width=5,
              withSpinner(DTOutput("datasets_tbl"))
          ),
          column(
            width=7,
            fluidRow(
              column(
                width=6, 
                uiOutput("dataset_summary"),
              ),
              column(
                width=6, 
                withSpinner(plotOutput("dataset_preview_plot", height="800px")),
              )
            ),
          ),
        )
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
