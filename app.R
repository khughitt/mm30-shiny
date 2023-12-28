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
library(plotly)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(survminer)
library(svglite)
library(tidyverse)
library(yaml)

options(spinner.color = "#00bc8c")
options(scipen = 2, digits = 3)
set.seed(1)

# table options
tableOpts <- list(pageLength = 15)

# TESTING (Dec23)
data_dir <- "../data"

# load config
cfg <- read_yaml("config/config-v7.0.yml")

# options for survival plot upper/lower expression cutoffs
# surv_expr_cutoffs <- cfg$surv_cutoff_opts
# names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")

# load mm30 gene scores
mm30_genes <- read_feather(file.path(data_dir, "mm30", "mm30_gene_scores.feather"))

# load metadata
gene_mdata <- read_tsv(file.path(data_dir, "metadata/genes.tsv"), col_types = cols()) %>%
  select(-chr_subband)

covariates <- read_yaml(file.path(data_dir, "metadata/covariates.yml"))

col_mdata <- readRDS(file.path(data_dir, "metadata/column-metadata.rds"))

# load disease stage -> sample id mapping
stage_sample_ids <- readRDS(file.path(data_dir, "metadata/disease_stage_sample_ids.rds"))

# pre-load most common tables; also helps by ensuring that genes available for each table are known
# at the time select elements are rendered
gene_rankings <- mm30_genes %>%
    left_join(gene_mdata, by = "symbol") %>%
    select(-num_datasets)

stage_rankings <- gene_rankings %>%
    select(Gene = symbol, Pval = disease_stage, Description = description,
           CHR = chr_region, `Cell Cycle` = cell_cycle_phase,
           Missing = num_missing_disease_stage) %>%
    filter(Missing <= cfg$max_missing$disease_stage) %>%
    arrange(Pval) %>%
    mutate(Rank = dense_rank(Pval)) %>%
    select(Rank, everything())

treatment_response_rankings <- gene_rankings %>%
    select(Gene = symbol, Pval = treatment_response, Description = description,
           CHR = chr_region, `Cell Cycle` = cell_cycle_phase,
           Missing = num_missing_treatment_response) %>%
    filter(Missing <= cfg$max_missing$treatment_response) %>%
    arrange(Pval) %>%
    mutate(Rank = dense_rank(Pval)) %>%
    select(Rank, everything())

# load mm30 individual GEO/MMRF gene expression datasets
geo_dir <- file.path(data_dir, "geo")

# slow; pre-filter genes and/or load indiv datasets on-the-fly?
#mm30_expr <- readRDS(file.path(data_dir, "expr/mm30_expr.rds"))

# load combined expr? (still need indiv expr datasets to access genes which are filtered out from
# the combined dataset, and for non-cpm transformed expr..)
combined_expr <- read_feather(file.path(data_dir, "expr/combined_expr_scaled.feather")) %>%
  column_to_rownames("symbol")

# load dataset + covariate metadata
covariate_mdata <- read_feather(file.path(data_dir, "metadata/covariates.feather"))

# individual dataset p-values
mm30_gene_pvals_indiv <- read_feather(file.path(data_dir,
                                                "fassoc/gene_association_pvals.feather"))


message("===========================================")
message("Initializing..")
message("===========================================")

#############################
#
# Shiny Server
#
#############################
server <- function(input, output, session) {
  #
  # Creates a violin plot for a given categorical response variable.
  #
  plot_categorical <- function(dat, dataset, covariate, feat_name, color_pal) {
    # drop any entries with missing values
    dat <- dat[!is.na(dat$response), ]

    # draw violin + jitter plot
    set.seed(1)

    ncat <- nlevels(dat$response)

    # if more factor levels exist than colors, expand palette
    if (ncat > length(color_pal)) {
      color_pal <- colorRampPalette(color_pal)(nlevels(dat$response))
    }

    plt <- ggplot(dat, aes(x = response, y = feature)) +
      geom_violin(aes(fill = response, color = response), alpha = 0.5, draw_quantiles = c(0.5)) +
      geom_boxplot(width = 0.1) +
      geom_jitter(aes(color = response), alpha = 0.8) +
      scale_fill_manual(values = color_pal) +
      scale_color_manual(values = color_pal, guide = "none") +
      theme_pubr(base_size = 16) +
      xlab(covariate) +
      ylab(sprintf("%s expression", feat_name)) +
      ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, covariate))

    # work-around to hide part of legend with plotly/ggplotly
    # https://github.com/plotly/plotly.R/issues/572#issuecomment-876540008
    plt %>%
      style(plt, showlegend = FALSE, traces = (ncat + 1):(ncat * 2 + 1))
  }

  # Creates a scatter plot for differential expression covariates
  plot_deseq <- function(dat, dataset, covariate, feat_name, color_pal) {
    # drop any entries with missing values
    #dat <- dat[!is.na(dat$response), ]

    # if more factor levels exist than colors, expand palette
    if (nlevels(dat$response) > length(color_pal)) {
      color_pal <- colorRampPalette(color_pal)(nlevels(dat$response))
    }

    var1 <- colnames(dat)[2]

    # single variable models
    if (ncol(dat) == 2) {
      ggplot(dat) +
        geom_col(aes(x = .data[[var1]], y = feature, fill = .data[[var1]], color = .data[[var1]])) +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, var1)) +
        theme_pubr(base_size = 16) +
        xlab(covariate) +
        ylab(sprintf("%s expression", feat_name))
    }
  }

    #
    # datasets
    #
    gse9782_expr <- reactive({
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

  #
  # html output
  #

  #
  # tables
  #

  output$disease_stage_tbl <- renderDataTable({
    #req(gene_df())
    DT::datatable(stage_rankings, style = "bootstrap", escape = FALSE,  rownames = FALSE,
                  selection = "single", options = tableOpts) %>%
      formatSignif(columns = c("Pval"), digits = 3)
  })

  output$treatment_response_tbl <- renderDataTable({
    #req(gene_df())
    DT::datatable(treatment_response_rankings,
                  style = "bootstrap", escape = FALSE,  rownames = FALSE,
                  selection = "single", options = tableOpts) %>%
      formatSignif(columns = c("Pval"), digits = 3)
  })

  #
  # plots
  #
  output$stage_plot <- renderPlot({
    set.seed(1)

    # create dataframe with disease stage gene expr estimates
    gene_expr <- combined_expr[input$gene_stg, , drop = FALSE]

    stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

    stage_expr_lst <- list(
      "Healthy" = c(),
      "MGUS" = c(),
      "SMM" = c(),
      "MM" = c(),
      "RRMM" = c()
      )

    for (stage in stages) {
      stage_ids <- stage_sample_ids[[stage]]
      stage_expr_vals <- as.numeric(gene_expr[, stage_ids])
      stage_expr_vals <- stage_expr_vals[!is.na(stage_expr_vals)]
      stage_expr_lst[[stage]] <- stage_expr_vals
    }

    # create long df
    stage_df <- stack(stage_expr_lst) %>%
      select(stage = ind, expr = values)

    # draw violin + jitter plot
    ggplot(stage_df, aes(x = stage, y = expr)) +
      geom_violin(aes(fill = stage, color = stage), alpha = 0.5, draw_quantiles = c(0.5)) +
      geom_jitter(aes(color = stage), alpha = 0.8) +
      scale_fill_manual(values = cfg$colors) +
      scale_color_manual(values = cfg$colors) +
      theme_pubr(base_size = 16)
      #ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, covariate)) +
      #xlab(covariate) +
      #ylab(sprintf("%s expression", feat_name))
  })

  output$gse9782_plot <- renderPlotly({
    gene_expr <- gse9782_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse9782_mdat() %>%
      select(geo_accession, response = treatment_response)

    df$feature <- gene_expr

    plot_categorical(df, "GSE9782", "Treatment Response (VTD)", input$gene_trmt, cfg$colors)
  })

  output$gse68871_plot <- renderPlotly({
    gene_expr <- gse68871_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse68871_mdat() %>%
      select(geo_accession, response = treatment_response)

    df$feature <- gene_expr

    plot_categorical(df, "GSE68871", "Treatment Response (VTD)", input$gene_trmt, cfg$colors)
  })

  output$gse39754_plot <- renderPlotly({
    gene_expr <- gse39754_expr() %>%
      filter(symbol == input$gene_trmt) %>%
      select(-symbol) %>%
      as.numeric()

    df <- gse39754_mdat() %>%
      select(geo_accession, response = treatment_response)

    df$feature <- gene_expr

    plot_categorical(df, "GSE39754", "Treatment Response (VAD + ACST)", input$gene_trmt, cfg$colors)
  })

  output$treatment_response_plots <- renderUI({
    plts <- list()

    if (input$gene_trmt %in% gse9782_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("gse9782_plot")))
    }

    if (input$gene_trmt %in% gse68871_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("gse68871_plot")))
    }

    if (input$gene_trmt %in% gse39754_expr()$symbol) {
      plts <- c(plts, list(plotlyOutput("gse39754_plot")))
    }

    tagList(plts)
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
      id = "tab",
      theme = shinytheme("flatly"),
      title = textOutput("page_title"),
      windowTitle = "MM30",

      tabPanel(
        "Disease Stage",
        fluidRow(
          column(
              width = 4,
              withSpinner(dataTableOutput("disease_stage_tbl"))
          ),
          column(
              width = 4,
              selectizeInput("gene_stg", "Gene:", choices = stage_rankings$Gene),
              helpText("Gene to visualize disease stage data for."),
              hr(),
              withSpinner(plotOutput("stage_plot", height = "740px"))
          )
        )
      ),
      tabPanel(
        "Survival",
        ""
      ),
      tabPanel(
        "Treatment",
        fluidRow(
          column(
              width = 4,
              withSpinner(dataTableOutput("treatment_response_tbl"))
          ),
          column(
              width = 4,
              fluidRow(
                column(
                  width = 4,
                  selectizeInput("gene_trmt", "Gene:",
                                 choices = treatment_response_rankings$Gene),
                  helpText("Gene to visualize treatment response data for.")
                ),
                column(
                 width = 2,
                 switchInput("show_indiv_treatment_plots", value = TRUE),
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
        ""
      ),
      tabPanel(
        "Co-Expression",
        ""
      ),
      tabPanel(
        "NLP",
        ""
      ),
      tabPanel(
        "Datasets",
        ""
      )
    )
  )
}

#
# Creates a kaplan-meyer survival plot for specified upper/lower expression quantiles.
#
plot_survival <- function(dat, dataset, covariate, feat_name, expr_cutoffs, color_pal) {
  # divide expression into quantiles
  cutoff <- as.numeric(expr_cutoffs)

  expr_quantiles <- quantile(dat$feature, c(cutoff / 100, 1 - (cutoff / 100)))

  dat$Expression <- ""
  dat$Expression[dat$feature <= expr_quantiles[1]] <- paste0("Lower ", names(expr_quantiles)[1])
  dat$Expression[dat$feature >= expr_quantiles[2]] <- paste0("Upper ", names(expr_quantiles)[1])

  # determine units to use
  if (dataset %in% c("MMRF", "GSE7039", "GSE57317", "GSE9782")) {
    time_units <- "days"
  } else if (dataset %in% c("GSE24080")) {
    time_units <- "weeks"
  } else if (dataset %in% c("GSE19784")) {
    time_units <- "months"
  } else {
    time_units <- "?"
  }

  # drop all data except for upper and lower quantiles
  dat <- dat %>%
    filter(Expression != "")

  num_samples <- nrow(dat)

  dat$Expression <- factor(dat$Expression)

  cov_label <- str_to_title(gsub("_", " ", covariate))
  plt_title <- sprintf("%s: %s vs. %s (n = %d)", dataset, cov_label, feat_name, num_samples)

  # perform fit on binarized data
  fit <- survival::survfit(survival::Surv(time, event) ~ Expression, data = dat)

  # display a kaplan meier plot for result
  survminer::ggsurvplot(fit, data = dat, ggtheme = theme_pubr(base_size = 16), palette = color_pal,
                        title = plt_title, xlab = sprintf("Time (%s)", time_units),
                        legend = "bottom", legend.title = "Legend")

}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
