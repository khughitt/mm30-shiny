#
# MM30 Shiny App
#
library(annotables)
library(arrow)
library(DT)
library(gridExtra)
library(ggpubr)
library(heatmaply)
library(Hmisc)
library(logger)
library(plotly)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(survminer)
library(svglite)
library(tidyverse)
library(yaml)
source("R/data.R")
source("R/plotting.R")

#############################
#
# Setup
#
#############################
options(spinner.color="#00bc8c")
options(scipen=2, digits=3)
set.seed(1)

log_threshold(DEBUG)
log_info("Initializing MM30 shiny app..")

# DT options
tableOpts <- list(pageLength=15, scrollX=TRUE)
selectOpts <- list(target="row", mode="single", selected=1)

# load config
cfg <- read_yaml("config/config-v7.3.yml")

results_dir <- file.path(cfg$data_dir, "results")
fassoc_dir <- file.path(cfg$data_dir, "fassoc")
mmrf_dir <- file.path(cfg$data_dir, "mmrf", "IA22")
geo_dir <- file.path(cfg$data_dir, "geo")

# ordered vector of disease stages
# disease_stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

# options for survival plot upper/lower expression cutoffs
# surv_expr_cutoffs <- cfg$surv_cutoff_opts
# names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")
surv_expr_cutoffs <- c(25, 75)

#############################
#
# Load data
#
#############################

log_info("Loading data..")

# load metadata
gene_mdata <- read_feather(file.path(results_dir, "metadata/genes.feather")) %>%
  select(-chr_subband)

sample_mdata <- read_feather(file.path(results_dir, "metadata/samples.feather"))
covariate_mdata <- read_feather(file.path(results_dir, "metadata/covariates.feather"))
indiv_mdata <- get_sample_metadata_list(mmrf_dir, geo_dir)

# load overall survival data
gene_scores_all <- read_feather(file.path(results_dir, "scores/metap/all/gene.feather")) %>%
  left_join(gene_mdata, by='symbol')

gene_scores_all$Rank <- 1:nrow(gene_scores_all)

gene_scores_surv_os <- read_feather(file.path(results_dir, "scores/combined/survival_os/gene.feather")) %>%
    left_join(gene_mdata, by='symbol')

surv_os_num_datasets <- max(gene_scores_surv_os$num_present)

gene_scores_surv_os <- gene_scores_surv_os %>%
    select(Gene=symbol, sumz_wt_pval, `P-value\n(metafor)`=metafor_pval, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase, Missing=num_missing) %>%
    filter(Missing <= cfg$max_missing$surv_os) %>%
    #rename(!!sprintf("Missing (N/%d)", surv_os_num_datasets) := Missing) %>%
    select(-Missing) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

surv_os_pvals <- read_feather(file.path(results_dir, "associations/survival_os/gene/pvals.feather"))
surv_os_errors <- read_feather(file.path(results_dir, "associations/survival_os/gene/errors.feather"))
surv_os_effects <- read_feather(file.path(results_dir, "associations/survival_os/gene/effects.feather"))

# load progress-free survival data
gene_scores_surv_pfs <- read_feather(file.path(results_dir, "scores/combined/survival_pfs/gene.feather")) %>%
    left_join(gene_mdata, by='symbol')

surv_pfs_num_datasets <- max(gene_scores_surv_pfs$num_present)

gene_scores_surv_pfs <- gene_scores_surv_pfs %>%
    select(Gene=symbol, sumz_wt_pval, `P-value\n(metafor)`=metafor_pval, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase, Missing=num_missing) %>%
    filter(Missing <= cfg$max_missing$surv_pfs) %>%
    #rename(!!sprintf("Missing (N/%d)", surv_pfs_num_datasets) := Missing) %>%
    select(-Missing) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

# load combined expr
expr_dat <- read_feather(file.path(results_dir, "expr/gene/expr_scaled.feather")) %>%
  column_to_rownames("symbol")

# load TCGA survival data
tcga <- read_feather("../scripts/data/hpa/hpa-cancer.feather")

# load HMCL gene expr data
expr_hmcl <- read_feather("../data/7.3/expr/expr_hmcl.feather")
hmcl_cells <- colnames(expr_hmcl)[-1]

log_info("Finished loading data..")

#############################
#
# Shiny Server
#
#############################
server <- function(input, output, session) {
  log_info("shiny::server()")

  ##########
  #
  # Overall survival
  #
  ###########
  surv_os_datasets <- covariate_mdata %>% 
    filter(phenotype=='overall_survival') %>% 
    pull(dataset) %>%
    sort(TRUE)

  surv_os_gene_selected <- reactive({
    gene_scores_surv_os$Gene[input$surv_os_tbl_rows_selected]
  })

  surv_os_gene_tcga_tbl <- reactive({
    req(input$surv_os_tbl_rows_selected)
    gene <- surv_os_gene_selected()

    tcga %>%
      filter(symbol == gene) %>%
      select(-symbol, -ensgene)
  })

  surv_os_gene_details_tbl <- reactive({
    req(input$surv_os_tbl_rows_selected)
    gene <- surv_os_gene_selected()

    df <- bind_rows(
      surv_os_pvals %>%
        filter(symbol == gene),
      surv_os_effects %>%
        filter(symbol == gene),
      surv_os_errors %>%
        filter(symbol == gene)
    )

    df %>%
      select(-symbol) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("P-value", "Hazard Ratio", "Std. Error")) %>%
      rownames_to_column("Dataset") %>%
      mutate(Dataset=str_replace(Dataset, '_overall_survival', ''))
  })

  output$surv_os_gene_summary_html <- renderUI({
      req(input$surv_os_tbl_rows_selected)
      gene <- surv_os_gene_selected()

      all_rank <- gene_scores_all %>% 
        filter(symbol == gene) %>%
        pull(Rank)

      surv_os_rank <- gene_scores_surv_os %>% 
        filter(Gene == gene) %>%
        pull(Rank)

      surv_pfs_rank <- gene_scores_surv_pfs %>% 
        filter(Gene == gene) %>%
        pull(Rank)

      tagList(
        tags$h2(gene),
        br(),
        tags$h4("MM30 Rankings:"),
        br(),
        HTML(sprintf("All: <b>%d</b>",  all_rank)), 
        br(),
        HTML(sprintf("Overall Survival (OS): <b>%d</b>", surv_os_rank)), 
        br(),
        HTML(sprintf("Progression-free Survival (PFS): <b>%d</b>", surv_pfs_rank))
      )
  })

  output$surv_os_details_tbl <- renderTable({
    surv_os_gene_details_tbl()
  }, digits=-3)

  output$surv_os_gene_tcga_tbl <- renderTable({
    surv_os_gene_tcga_tbl()
  }, digits=4)

  output$surv_os_tbl <- renderDT({
    DT::datatable(gene_scores_surv_os, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$surv_os_plot <- renderPlot({
    req(input$surv_os_tbl_rows_selected)

    plts <- list()

    # lapply(surv_os_datasets, function(acc) {
    for (acc in surv_os_datasets) {

      df <- indiv_mdata[[acc]] %>%
        select(1, time=os_time, event=os_censor)

      sample_ids <- df %>%
        pull(1)

      gene_name <- surv_os_gene_selected()
      gene_expr <- as.numeric(expr_dat[gene_name, sample_ids])

      df$feature <- gene_expr 

      time_units <- "Days"

      if (acc %in% c("GSE106218", "GSE19784", "GSE24080")) {
        time_units <- "Months"
      }

      if (sum(!is.na(gene_expr)) > 0) {
        plts[[acc]] <- plot_survival(df, acc, time_units, surv_expr_cutoffs, cfg$colors)
      }

    # })
    }

    arrange_ggsurvplots(plts, nrow=ceiling(length(plts) / 3), ncol=3,
                        title=sprintf("Overall survival (%s)", gene_name))
  })

  ##########
  #
  # Cell lines (HMCL)
  #
  ###########
  output$hmcl_expr_tbl <- renderDT({
    DT::datatable(expr_hmcl, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=hmcl_cells, digits=3)
  })

  hmcl_gene_selected <- reactive({
    expr_hmcl$symbol[input$hmcl_expr_tbl_rows_selected]
  })

  output$hmcl_hist_plot <- renderPlot({
    req(input$hmcl_expr_tbl_rows_selected)

    gene <- hmcl_gene_selected()

    x <- expr_hmcl %>%
      filter(symbol == gene) %>%
      select(-symbol) %>%
      as.numeric()

    df <- data.frame(x=x)

    ggplot(df, aes(x)) +
      geom_histogram() +
      ggtitle(sprintf("HMCL Expression (%s)", gene)) +
      xlab("Expression (TPM)") +
      ylab("Count")
  })

  # bookmarking support
  observe({
    reactiveValuesToList(input)
    session$doBookmark()
  })
  onBookmarked(updateQueryString)

  setBookmarkExclude(c(
    "disease_stage_tbl_rows_selected",
    "disease_stage_tbl_columns_selected",
    "disease_stage_tbl_cells_selected",
    "disease_stage_tbl_rows_current",
    "disease_stage_tbl_rows_all",
    "disease_stage_tbl_state",
    "disease_stage_tbl_search",
    "disease_stage_tbl_cell_clicked",
    "treatment_response_tbl_rows_selected",
    "treatment_response_tbl_columns_selected",
    "treatment_response_tbl_cells_selected",
    "treatment_response_tbl_rows_current",
    "treatment_response_tbl_rows_all",
    "treatment_response_tbl_state",
    "treatment_response_tbl_search",
    "treatment_response_tbl_cell_clicked"
  ))
}

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
      tabPanel(
        "OS",
        fluidRow(
          column(
            width=6,
            withSpinner(DTOutput("surv_os_tbl"))
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
            tableOutput("surv_os_details_tbl"),
            tagList(
              tags$hr(),
              tags$h3("TCGA survival associations:")
            ),
            tableOutput("surv_os_gene_tcga_tbl")
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
            #withSpinner(plotOutput("surv_os_plot", height="800px"))
          ),
        ),
      ),
      tabPanel(
        "[Co-Expression]",
        ""
      ),
      tabPanel(
        "[NLP?]",
        ""
      ),
      tabPanel(
        "[Datasets]",
        ""
      )
    )
  )
}

shinyApp(ui=ui, server=server, enableBookmarking="url")
