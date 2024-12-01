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
source("R/plotting.R")

options(spinner.color="#00bc8c")
options(scipen=2, digits=3)
set.seed(1)

log_threshold(DEBUG)
log_info("Initializing MM30 shiny app..")

tableOpts <- list(pageLength=15)

data_dir <- "../data"

# load config
cfg <- read_yaml("config/config-v7.0.yml")

# number of datasets used for each subset
num_datasets <- list(
  "stage"=17,
  "surv_os"=6,
  "surv_pfs"=4,
  "treatment"=4
)

# ordered vector of disease stages
# disease_stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

# options for survival plot upper/lower expression cutoffs
# surv_expr_cutoffs <- cfg$surv_cutoff_opts
# names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")
surv_expr_cutoffs <- c(25, 75)

# load mm30 gene scores
log_info("Loading mm30_gene_scores.feather")
mm30_genes <- read_feather(file.path(data_dir, "mm30", "mm30_gene_scores.feather"))

# load metadata
gene_mdata <- read_tsv(file.path(data_dir, "metadata/genes.tsv"), show_col_types=FALSE) %>%
  select(-chr_subband)

covariates <- read_yaml(file.path(data_dir, "metadata/covariates.yml"))

col_mdata <- readRDS(file.path(data_dir, "metadata/column-metadata.rds"))

# load disease stage -> sample id mapping
stage_sample_ids <- readRDS(file.path(data_dir, "metadata/disease_stage_sample_ids.rds"))

# load HMCL CPM count data
hmcl_expr <- read_feather(file.path(data_dir, "hmcl", "rnaseq_symbols.feather"))

# pre-load most common tables; also helps by ensuring that genes available for each table are known
# at the time select elements are rendered
log_info("Creating gene ranking dataframes")

gene_rankings <- mm30_genes %>%
    left_join(gene_mdata, by="symbol") %>%
    select(-num_datasets)

stage_rankings <- gene_rankings %>%
    select(Gene=symbol, Pval=disease_stage, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase,
           DGIdb=dgidb_categories,
           Missing=disease_stage_num_missing) %>%
    filter(Missing <= cfg$max_missing$disease_stage) %>%
    rename(!!sprintf("Missing (N / %d)", num_datasets$stage) := Missing) %>%
    arrange(Pval) %>%
    mutate(Rank=dense_rank(Pval)) %>%
    select(Rank, everything())

# highlight genes located on ch1
ch1 <- substr(stage_rankings$CHR, 1, 2) %in% c("1p", "1q")
stage_rankings$chr_color <- ifelse(ch1, "#77BDF3", "")

survival_os_rankings <- gene_rankings %>%
    select(Gene=symbol, Pval=survival_os, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase,
           Missing=survival_os_num_missing) %>%
    filter(Missing <= cfg$max_missing$surv_os) %>%
    rename(!!sprintf("Missing (N/%d)", num_datasets$surv_os) := Missing) %>%
    arrange(Pval) %>%
    mutate(Rank=dense_rank(Pval)) %>%
    select(Rank, everything())

survival_pfs_rankings <- gene_rankings %>%
    select(Gene=symbol, Pval=survival_pfs, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase,
           Missing=survival_pfs_num_missing) %>%
    filter(Missing <= cfg$max_missing$surv_pfs) %>%
    rename(!!sprintf("Missing (N/%d)", num_datasets$surv_pfs) := Missing) %>%
    arrange(Pval) %>%
    mutate(Rank=dense_rank(Pval)) %>%
    select(Rank, everything())

treatment_response_rankings <- gene_rankings %>%
    select(Gene=symbol, Pval=treatment_response, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase,
           Missing=treatment_response_num_missing) %>%
    filter(Missing <= cfg$max_missing$treatment_response) %>%
    rename(!!sprintf("Missing (N/%d)", num_datasets$treatment) := Missing) %>%
    arrange(Pval) %>%
    mutate(Rank=dense_rank(Pval)) %>%
    select(Rank, everything())

# load mm30 individual GEO/MMRF gene expression datasets
geo_dir <- file.path(data_dir, "geo")
mmrf_dir <- file.path(data_dir, "mmrf")

# slow; pre-filter genes and/or load indiv datasets on-the-fly?
#mm30_expr <- readRDS(file.path(data_dir, "expr/mm30_expr.rds"))

# load combined expr? (still need indiv expr datasets to access genes which are filtered out from
# the combined dataset, and for non-cpm transformed expr..)
log_info("loading combined_expr_scaled.feather")

combined_expr <- read_feather(file.path(data_dir, "expr/combined_expr_scaled.feather")) %>%
  column_to_rownames("symbol")

# load dataset + covariate metadata
covariate_mdata <- read_feather(file.path(data_dir, "metadata/covariates.feather"))

# individual dataset p-values
log_info("loading gene_association_pvals.feather")
mm30_gene_pvals_indiv <- read_feather(file.path(data_dir,
                                                "fassoc/gene_association_pvals.feather"))

#############################
#
# Shiny Server
#
#############################
server <- function(input, output, session) {
  log_info("shiny::server()")

  ################################################################################
  #
  # datasets
  #
  ################################################################################

  # GEO data
  gse106218_expr <- reactive({
    log_info("loading GSE106218/data.feather")
    read_feather(file.path(geo_dir, "GSE106218", "data.feather"))
  })
  gse106218_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE106218", "column-metadata.feather"))
    df$iss_stage <- ordered(df$iss_stage, c("I", "II", "III"))
    df
  })

  gse117846_expr <- reactive({
    log_info("loading GSE117846/data.feather")
    read_feather(file.path(geo_dir, "GSE117846", "data.feather"))
  })
  gse117846_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE117846", "column-metadata.feather"))
    df$disease_stage <- ordered(df$disease_stage, c("SMM", "MM"))
    df
  })

  gse118900_expr <- reactive({
    log_info("loading GSE118900/data.feather")
    read_feather(file.path(geo_dir, "GSE118900", "data.feather"))
  })
  gse118900_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE118900", "column-metadata.feather"))
    df$disease_stage <- ordered(df$disease_stage, c("MGUS", "SMM", "MM", "RRMM"))
    df
  })

  gse128251_expr <- reactive({
    log_info("loading GSE128251/data.feather")
    read_feather(file.path(geo_dir, "GSE128251", "data.feather"))
  })
  gse128251_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE128251", "column-metadata.feather"))
    df$treatment <- factor(df$treatment)
    df$replicate <- factor(df$replicate)
    df
  })

  gse134598_expr <- reactive({
    log_info("loading GSE134598/data.feather")
    read_feather(file.path(geo_dir, "GSE134598", "data.feather"))
  })
  gse134598_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE134598", "column-metadata.feather"))
    df$cell_line <- factor(df$cell_line)
    df$treatment <- factor(df$treatment)
    df$dose <- as.numeric(sub("nM", "", df$dose))
    df
  })

  gse13591_expr <- reactive({
    log_info("loading GSE13591/data.feather")
    read_feather(file.path(geo_dir, "GSE13591", "data.feather"))
  })
  gse13591_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE13591", "column-metadata.feather"))
    df$disease_stage <- ordered(df$disease_stage, c("Healthy", "MGUS", "MM"))
    df
  })

  gse144249_expr <- reactive({
    log_info("loading GSE144249/data.feather")
    read_feather(file.path(geo_dir, "GSE144249", "data.feather"))
  })
  gse144249_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE144249", "column-metadata.feather"))
    df$drug_resistance <- factor(df$drug_resistance)
    df$replicate <- factor(df$replicate)
    df
  })

  ########

  gse19784_expr <- reactive({
    log_info("loading GSE19784/data.feather")
    read_feather(file.path(geo_dir, "GSE19784", "data.feather"))
  })
  gse19784_mdat <- reactive({
    read_feather(file.path(geo_dir, "GSE19784", "column-metadata.feather"))
  })

  gse9782_expr <- reactive({
    log_info("loading GSE9782/data.feather")
    read_feather(file.path(geo_dir, "GSE9782", "data.feather"))
  })
  gse9782_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE9782", "column-metadata.feather"))
    df$treatment_response <- ordered(df$treatment_response, c("CR", "PR", "MR", "NC", "PD"))
    df
  })

  gse68871_expr <- reactive({
    read_feather(file.path(geo_dir, "GSE68871", "data.feather"))
  })
  gse68871_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE68871", "column-metadata.feather"))
    df$treatment_response <- ordered(df$treatment_response, c("CR", "nCR", "VGPR", "PR", "SD"))
    df
  })

  gse39754_expr <- reactive({
    read_feather(file.path(geo_dir, "GSE39754", "data.feather"))
  })
  gse39754_mdat <- reactive({
    df <- read_feather(file.path(geo_dir, "GSE39754", "column-metadata.feather"))
    df$treatment_response <- ordered(df$treatment_response, c("CR", "VGPR", "PR", "NR", "Prog"))
    df
  })

  # MMRF data
  mmrf_expr <- reactive({
    log_info("loading mmrf/data.feather")
    read_feather(file.path(mmrf_dir, "data.feather"))
  })
  mmrf_mdat <- reactive({
    df <- read_feather(file.path(mmrf_dir, "column-metadata.feather")) %>%
      select(public_id, iss_stage, ecog, mm_status, fresp, frespcd,
              trtshnm, pfscdy, censpfs, oscdy, censos)

    df$frespcd <- ordered(df$frespcd, c("sCR", "CR", "VGPR", "PR"))

    df$response_bor_len_dex <- df$frespcd
    df$response_bor_cyc_dex <- df$frespcd

    df$response_bor_len_dex[df$trtshnm != "Bor-Len-Dex"] <- NA
    df$response_bor_cyc_dex[df$trtshnm != "Bor-Cyc-Dex"] <- NA

    df
  })

  ################################################################################
  #
  # tables
  #
  ################################################################################
  output$hmcl_expr_tbl <- renderDataTable({
    dat <- hmcl_expr  %>%
      filter(symbol == input$gene_hmcl) %>%
      select(-symbol) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("cell_line") %>%
      rename(expr=V1) %>%
      arrange(-expr)

    DT::datatable(dat, style="bootstrap", options=tableOpts) %>%
      formatSignif(columns=c("expr"), digits=3)
  })

  output$disease_stage_tbl <- renderDataTable({
    tblOpts <- c(tableOpts, list(columnDefs=list(list(visible=FALSE, targets=c("chr_color")))))

    DT::datatable(stage_rankings,
                  style="bootstrap",
                  escape=FALSE,
                  rownames=FALSE,
                  selection="single",
                  options=tblOpts) %>%
      formatSignif(columns=c("Pval"), digits=3) %>%
      formatStyle(
            "CHR",
            valueColumns="chr_color",
            backgroundColor=JS("value")
          )

  })

  output$survival_tbl <- renderDataTable({
    if (input$surv_subset == "OS") {
      dat <- survival_os_rankings
    } else {
      dat <- survival_pfs_rankings
    }

    DT::datatable(dat, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection="single", options=tableOpts) %>%
      formatSignif(columns=c("Pval"), digits=3)
  })

  output$treatment_response_tbl <- renderDataTable({
    DT::datatable(treatment_response_rankings,
                  style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection="single", options=tableOpts) %>%
      formatSignif(columns=c("Pval"), digits=3)
  })

  ################################################################################
  #
  # disease stage plots
  #
  ################################################################################
  output$stage_plot <- renderPlot({
    set.seed(1)
    log_info("output$stage_plot")

    # create dataframe with disease stage gene expr estimates
    gene_expr <- combined_expr[input$gene_stg, , drop=FALSE]

    stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

    stage_expr_lst <- list(
      "Healthy"=c(),
      "MGUS"=c(),
      "SMM"=c(),
      "MM"=c(),
      "RRMM"=c()
      )

    for (stage in stages) {
      stage_ids <- stage_sample_ids[[stage]]
      stage_expr_vals <- as.numeric(gene_expr[, stage_ids])
      stage_expr_vals <- stage_expr_vals[!is.na(stage_expr_vals)]
      stage_expr_lst[[stage]] <- stage_expr_vals
    }

    # create long df
    stage_df <- stack(stage_expr_lst) %>%
      select(stage=ind, expr=values)

    # draw violin + jitter plot
    ggplot(stage_df, aes(x=stage, y=expr)) +
      geom_violin(aes(fill=stage, color=stage), alpha=0.5, draw_quantiles=c(0.5)) +
      geom_jitter(aes(color=stage), alpha=0.8) +
      scale_fill_manual(values=cfg$colors) +
      scale_color_manual(values=cfg$colors) +
      theme_pubr(base_size=16)
      #ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, covariate)) +
      #xlab(covariate) +
      #ylab(sprintf("%s expression", feat_name))
  })

  ################################################################################
  #
  # treatment response plots
  #
  ################################################################################
  output$gse9782_treatment_plot <- renderPlotly({
    log_info("output$gse9782_treatment_plot")

    gene_expr <- gse9782_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse9782_mdat() %>%
      select(geo_accession, response=treatment_response)

    df$feature <- gene_expr

    plot_categorical(df, "GSE9782", "Treatment Response (VTD)", input$gene_trmt, cfg$colors)
  })

  output$gse68871_treatment_plot <- renderPlotly({
    gene_expr <- gse68871_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse68871_mdat() %>%
      select(geo_accession, response=treatment_response)

    df$feature <- gene_expr

    plot_categorical(df, "GSE68871", "Treatment Response (VTD)", input$gene_trmt, cfg$colors)
  })

  output$gse39754_treatment_plot <- renderPlotly({
    gene_expr <- gse39754_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse39754_mdat() %>%
      select(geo_accession, response=treatment_response)

    df$feature <- gene_expr

    plot_categorical(df, "GSE39754", "Treatment Response (VAD + ACST)", input$gene_trmt, cfg$colors)
  })

  output$mmrf_bor_len_dex_plot <- renderPlotly({
    gene_expr <- mmrf_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- mmrf_mdat() %>%
      select(public_id, response=response_bor_len_dex)

    df$feature <- gene_expr

    plot_categorical(df, "MMRF", "Treatment Response (Bor-Len-Dex)", input$gene_trmt, cfg$colors)
  })

  output$mmrf_bor_cyc_dex_plot <- renderPlotly({
    gene_expr <- mmrf_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- mmrf_mdat() %>%
      select(public_id, response=response_bor_cyc_dex)

    df$feature <- gene_expr

    plot_categorical(df, "MMRF", "Treatment Response (Bor-Cyc-Dex)", input$gene_trmt, cfg$colors)
  })

  output$treatment_response_plots <- renderUI({
    log_info("output$treatment_response_plots")

    plts <- list()

    if (input$gene_trmt %in% mmrf_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("mmrf_bor_len_dex_plot")))
      plts <- c(plts, list(plotlyOutput("mmrf_bor_cyc_dex_plot")))
    }

    if (input$gene_trmt %in% gse9782_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("gse9782_treatment_plot", height="740px")))
    }

    if (input$gene_trmt %in% gse68871_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("gse68871_treatment_plot")))
    }

    if (input$gene_trmt %in% gse39754_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("gse39754_treatment_plot")))
    }

    tagList(plts)
  })

  ################################################################################
  #
  # survival plots
  #
  ################################################################################
  output$mmrf_os_plot <- renderPlot({
    # could also make a reactive "gene_expr_single_gene()" to allow reactivity to trigger?
    #req(input$gene_surv_os)

    log_info(sprintf("mmrf_os_plot(%s)", input$gene_surv_os))

    gene_expr <- mmrf_expr() %>%
      filter(symbol == input$gene_surv_os) %>%
      select(-symbol) %>%
      as.numeric()

    df <- mmrf_mdat() %>%
      select(public_id, time=oscdy, event=censos)

    df$feature <- gene_expr

    plot_survival(df, "MMRF", "Overall Survival", "Gene Expression", "days",
                  surv_expr_cutoffs, cfg$colors)
  })

  output$gse19784_os_plot <- renderPlot({
    gene_expr <- gse19784_expr() %>%
      filter(symbol == input$gene_surv_os) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse19784_mdat() %>%
      select(public_id, time=os_time, event=os_event)
    df$feature <- gene_expr

    plot_survival(df, "GSE19784", "Overall Survival", "Gene Expression", "Days",
                  surv_expr_cutoffs, cfg$colors)
  })

  output$gse19784_pfs_plot <- renderPlot({
    gene_expr <- gse19784_expr() %>%
      filter(symbol == input$gene_surv_pfs) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse19784_mdat() %>%
      select(public_id, time=pfs_time, event=pfs_event)
    df$feature <- gene_expr

    plot_survival(df, "GSE19784", "Progression Free Survival", "Gene Expression", "Days",
                  surv_expr_cutoffs, cfg$colors)
  })

  output$surv_os_plots <- renderUI({
    log_info("output$surv_os_plots")

    plts <- list()

    if (input$gene_surv_os %in% mmrf_expr()$symbol) {
      plts <- c(plts, list(plotOutput("mmrf_os_plot", height="740px")))
    }
  })

  # overall survival plots
  # MMRF
  # GSE19784
  # GSE24080
  # GSE7039
  # GSE57317
  # GSE9782

  # prog free survival plots
  # GSE9782
  # GSE19784
  # MMRF

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
        "Disease Stage",
        fluidRow(
          column(
              width=6,
              withSpinner(dataTableOutput("disease_stage_tbl"))
          ),
          column(
              width=6,
              selectizeInput("gene_stg", "Gene:", choices=stage_rankings$Gene),
              helpText("Gene to visualize disease stage data for."),
              hr(),
              withSpinner(plotOutput("stage_plot", height="740px"))
          )
        )
      ),
      tabPanel(
        "Survival",
        fluidRow(
                radioButtons("surv_subset", "OS/PFS?",
                             selected="OS",
                             inline=TRUE,
                             choiceNames=list("OS", "PFS"),
                             choiceValues=list("OS", "PFS"))
        ),
        fluidRow(
          column(
              width=6,
              withSpinner(dataTableOutput("survival_tbl"))
          ),
          column(
              width=6,
              fluidRow(
                column(
                  width=6,
                  selectizeInput("gene_surv_os", "Gene:",
                                 choices=survival_os_rankings$Gene),
                  helpText("Gene to visualize OS data for.")
                )
              ),
              hr(),
              withSpinner(htmlOutput("surv_os_plots"))
          )
        ),
      ),
      tabPanel(
        "Treatment",
        fluidRow(
          column(
              width=6,
              withSpinner(dataTableOutput("treatment_response_tbl"))
          ),
          column(
              width=6,
              fluidRow(
                column(
                  width=4,
                  selectizeInput("gene_trmt", "Gene:",
                                 choices=treatment_response_rankings$Gene),
                  helpText("Gene to visualize treatment response data for.")
                ),
                column(
                 width=2,
                 switchInput("show_indiv_treatment_plots", value=TRUE),
                 helpText("Show individual plots.")
                )
              ),
              hr(),
              withSpinner(htmlOutput("treatment_response_plots"))
          )
        )
      ),
      tabPanel(
        "Cell Lines",
        fluidRow(
            column(
              width=4,
              selectizeInput("gene_hmcl", "Gene:",
                              choices=hmcl_expr$symbol),
              helpText("Gene to show cell line gene expression data for.")
            ),
            column(
              width=6,
              withSpinner(dataTableOutput("hmcl_expr_tbl"))
            ),
        )
      ),
      tabPanel(
        "Co-Expression*",
        ""
      ),
      tabPanel(
        "CRISPR*",
        ""
      ),
      tabPanel(
        "NLP*",
        ""
      ),
      tabPanel(
        "Datasets*",
        ""
      )
    )
  )
}

shinyApp(ui=ui, server=server, enableBookmarking="url")
