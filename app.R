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

options(spinner.color="#00bc8c")
options(scipen=2, digits=3)
set.seed(1)

log_threshold(DEBUG)
log_info("Initializing MM30 shiny app..")

# DT options
tableOpts <- list(pageLength=15)
selectOpts <- list(target="row", mode="single", selected=1)

# load config
cfg <- read_yaml("config/config-v7.2.yml")

results_dir <- file.path(cfg$data_dir, "results")
fassoc_dir <- file.path(cfg$data_dir, "fassoc")
mmrf_dir <- file.path(cfg$data_dir, "mmrf", "IA22")
geo_dir <- file.path(cfg$data_dir, "geo")

# number of datasets used for each subset
# subset_sizes <- list(
#   "stage"=17,
#   "surv_os"=6,
#   "surv_pfs"=4,
#   "treatment"=4
# )

# ordered vector of disease stages
# disease_stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

# options for survival plot upper/lower expression cutoffs
# surv_expr_cutoffs <- cfg$surv_cutoff_opts
# names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")
surv_expr_cutoffs <- c(25, 75)

# load gene metadata
gene_mdata <- read_feather(file.path(results_dir, "metadata/genes.feather")) %>%
  select(-chr_subband)

#
# 1. Survival Rankings (OS / genes)
#
gene_scores_all <- read_feather(file.path(results_dir, "scores/metap/all/gene.feather")) %>%
  left_join(gene_mdata, by='symbol')

surv_os_gene_scores <- read_feather(file.path(results_dir, "scores/combined/survival_os/gene.feather")) %>%
    left_join(gene_mdata, by='symbol')

surv_os_num_datasets <- max(surv_os_gene_scores$num_present)

surv_os_gene_scores <- surv_os_gene_scores %>%
    select(Gene=symbol, `Metap P-value`=sumz_wt_pval, `Metafor P-value`=metafor_pval, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase, Missing=num_missing) %>%
    filter(Missing <= cfg$max_missing$surv_os) %>%
    rename(!!sprintf("Missing (N/%d)", surv_os_num_datasets) := Missing) %>%
    arrange(`Metap P-value`) %>%
    mutate(Rank=dense_rank(`Metap P-value`)) %>%
    select(Rank, everything())

#
# . Survival Rankings (PFS / genes)
#
surv_pfs_genes <- read_feather(file.path(results_dir, "scores/combined/survival_pfs/gene.feather")) %>%
    left_join(gene_mdata, by='symbol')

surv_pfs_num_datasets <- max(surv_pfs_genes$num_present)

surv_pfs_genes <- surv_pfs_genes %>%
    select(Gene=symbol, `Metap P-value`=sumz_wt_pval, `Metafor P-value`=metafor_pval, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase, Missing=num_missing) %>%
    filter(Missing <= cfg$max_missing$surv_pfs) %>%
    rename(!!sprintf("Missing (N/%d)", surv_pfs_num_datasets) := Missing) %>%
    arrange(`Metap P-value`) %>%
    mutate(Rank=dense_rank(`Metap P-value`)) %>%
    select(Rank, everything())


# load combined expr? (still need indiv expr datasets to access genes which are filtered out from
# the combined dataset, and for non-cpm transformed expr..)
log_info("loading expr_scaled.feather")

expr_dat <- read_feather(file.path(results_dir, "expr/gene/expr_scaled.feather")) %>%
  column_to_rownames("symbol")

# load sample metadata
sample_mdata <- read_feather(file.path(results_dir, "metadata/samples.feather"))

# load dataset + covariate metadata
covariate_mdata <- read_feather(file.path(results_dir, "metadata/covariates.feather"))

# individual dataset p-values
# log_info("loading gene_association_pvals.feather")
# mm30_gene_pvals_indiv <- read_feather(file.path(fassoc_dir,
#                                                 "merged/gene_association_pvals.feather"))

# 
# load overall survival data
#
surv_os_pvals <- read_feather(file.path(results_dir, "associations/survival_os/gene/pvals.feather"))
surv_os_errors <- read_feather(file.path(results_dir, "associations/survival_os/gene/errors.feather"))
surv_os_effects <- read_feather(file.path(results_dir, "associations/survival_os/gene/effects.feather"))

# load individual dataset sample metadata
indiv_mdata <- get_sample_metadata_list(mmrf_dir, geo_dir)

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
  surv_os_gene_selected <- reactive({
    surv_os_gene_scores$Gene[input$surv_os_tbl_rows_selected]
  })

  # surv_os_gene_summary <- reactive({
  # })

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
  ################################################################################
  #
  # datasets
  #
  ################################################################################

  # MMRF data
  # mmrf_expr <- reactive({
  #   log_info("loading mmrf/data.feather")
  #   read_feather(file.path(mmrf_dir, "data.feather"))
  # })

  ################################################################################
  #
  # tables
  #
  ################################################################################
  # output$hmcl_expr_tbl <- renderDT({
  #   dat <- hmcl_expr  %>%
  #     filter(symbol == input$gene_hmcl) %>%
  #     select(-symbol) %>%
  #     t() %>%
  #     as.data.frame() %>%
  #     rownames_to_column("cell_line") %>%
  #     rename(expr=V1) %>%
  #     arrange(-expr)
  #
  #   DT::datatable(dat, style="bootstrap", options=tableOpts) %>%
  #     formatSignif(columns=c("expr"), digits=3)
  # })
  #
  # output$disease_stage_tbl <- renderDT({
  #   tblOpts <- c(tableOpts, list(columnDefs=list(list(visible=FALSE, targets=c("chr_color")))))
  #
  #   DT::datatable(stage_rankings,
  #                 style="bootstrap",
  #                 escape=FALSE,
  #                 rownames=FALSE,
  #                 selection="single",
  #                 options=tblOpts) %>%
  #     formatSignif(columns=c("Metap P-value"), digits=3) %>%
  #     formatStyle(
  #           "CHR",
  #           valueColumns="chr_color",
  #           backgroundColor=JS("value")
  #         )
  #
  # })

  output$surv_os_details_tbl <- renderTable({
    surv_os_gene_details_tbl()
      # formatSignif(columns=c("Metap P-value", "Metafor P-value"), digits=3)
  })

  output$surv_os_tbl <- renderDT({
    DT::datatable(surv_os_gene_scores, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("Metap P-value", "Metafor P-value"), digits=3)
  })

  # output$treatment_response_tbl <- renderDT({
  #   DT::datatable(treatment_response_rankings,
  #                 style="bootstrap", escape=FALSE,  rownames=FALSE,
  #                 selection="single", options=tableOpts) %>%
  #     formatSignif(columns=c("Metap P-value"), digits=3)
  # })
  
  ################################################################################
  #
  # disease stage plots
  #
  ################################################################################
  # output$stage_plot <- renderPlot({
  #   set.seed(1)
  #   log_info("output$stage_plot")
  #
  #   # create dataframe with disease stage gene expr estimates
  #   gene_expr <- expr_dat[input$gene_stg, , drop=FALSE]
  #
  #   stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")
  #
  #   stage_expr_lst <- list(
  #     "Healthy"=c(),
  #     "MGUS"=c(),
  #     "SMM"=c(),
  #     "MM"=c(),
  #     "RRMM"=c()
  #     )
  #
  #   for (stage in stages) {
  #     stage_ids <- stage_sample_ids[[stage]]
  #     stage_expr_vals <- as.numeric(gene_expr[, stage_ids])
  #     stage_expr_vals <- stage_expr_vals[!is.na(stage_expr_vals)]
  #     stage_expr_lst[[stage]] <- stage_expr_vals
  #   }
  #
  #   # create long df
  #   stage_df <- stack(stage_expr_lst) %>%
  #     select(stage=ind, expr=values)
  #
  #   # draw violin + jitter plot
  #   ggplot(stage_df, aes(x=stage, y=expr)) +
  #     geom_violin(aes(fill=stage, color=stage), alpha=0.5, draw_quantiles=c(0.5)) +
  #     geom_jitter(aes(color=stage), alpha=0.8) +
  #     scale_fill_manual(values=cfg$colors) +
  #     scale_color_manual(values=cfg$colors) +
  #     theme_pubr(base_size=16)
  #     #ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, covariate)) +
  #     #xlab(covariate) +
  #     #ylab(sprintf("%s expression", feat_name))
  # })

  ################################################################################
  #
  # treatment response plots
  #
  ################################################################################
  # output$gse9782_treatment_plot <- renderPlotly({
  #   log_info("output$gse9782_treatment_plot")
  #
  #   gene_expr <- gse9782_expr() %>%
  #     filter(symbol == input$gene_trmt) %>%
  #     select(-symbol) %>%
  #     as.numeric()
  #
  #   df <- gse9782_mdat() %>%
  #     select(geo_accession, response=treatment_response)
  #
  #   df$feature <- gene_expr
  #
  #   plot_categorical(df, "GSE9782", "Treatment Response (VTD)", input$gene_trmt, cfg$colors)
  # })
  #
  # output$gse68871_treatment_plot <- renderPlotly({
  #   gene_expr <- gse68871_expr() %>%
  #     filter(symbol == input$gene_trmt) %>%
  #     select(-symbol) %>%
  #     as.numeric()
  #
  #   df <- gse68871_mdat() %>%
  #     select(geo_accession, response=treatment_response)
  #
  #   df$feature <- gene_expr
  #
  #   plot_categorical(df, "GSE68871", "Treatment Response (VTD)", input$gene_trmt, cfg$colors)
  # })
  #
  # output$gse39754_treatment_plot <- renderPlotly({
  #   gene_expr <- gse39754_expr() %>%
  #     filter(symbol == input$gene_trmt) %>%
  #     select(-symbol) %>%
  #     as.numeric()
  #
  #   df <- gse39754_mdat() %>%
  #     select(geo_accession, response=treatment_response)
  #
  #   df$feature <- gene_expr
  #
  #   plot_categorical(df, "GSE39754", "Treatment Response (VAD + ACST)", input$gene_trmt, cfg$colors)
  # })
  #
  # output$mmrf_bor_len_dex_plot <- renderPlotly({
  #   gene_expr <- mmrf_expr() %>%
  #     filter(symbol == input$gene_trmt) %>%
  #     select(-symbol) %>%
  #     as.numeric()
  #
  #   df <- mmrf_mdat() %>%
  #     select(public_id, response=response_bor_len_dex)
  #
  #   df$feature <- gene_expr
  #
  #   plot_categorical(df, "MMRF", "Treatment Response (Bor-Len-Dex)", input$gene_trmt, cfg$colors)
  # })
  #
  # output$mmrf_bor_cyc_dex_plot <- renderPlotly({
  #   gene_expr <- mmrf_expr() %>%
  #     filter(symbol == input$gene_trmt) %>%
  #     select(-symbol) %>%
  #     as.numeric()
  #
  #   df <- mmrf_mdat() %>%
  #     select(public_id, response=response_bor_cyc_dex)
  #
  #   df$feature <- gene_expr
  #
  #   plot_categorical(df, "MMRF", "Treatment Response (Bor-Cyc-Dex)", input$gene_trmt, cfg$colors)
  # })
  #
  # output$treatment_response_plots <- renderUI({
  #   log_info("output$treatment_response_plots")
  #
  #   plts <- list()
  #
  #   if (input$gene_trmt %in% mmrf_expr()$symbol) {
  #     plts <- c(plts, list(plotlyOutput("mmrf_bor_len_dex_plot")))
  #     plts <- c(plts, list(plotlyOutput("mmrf_bor_cyc_dex_plot")))
  #   }
  #
  #   if (input$gene_trmt %in% gse9782_expr()$symbol) {
  #     plts <- c(plts, list(plotlyOutput("gse9782_treatment_plot", height="740px")))
  #   }
  #
  #   if (input$gene_trmt %in% gse68871_expr()$symbol) {
  #     plts <- c(plts, list(plotlyOutput("gse68871_treatment_plot")))
  #   }
  #
  #   if (input$gene_trmt %in% gse39754_expr()$symbol) {
  #     plts <- c(plts, list(plotlyOutput("gse39754_treatment_plot")))
  #   }
  #
  #   tagList(plts)
  # })
  #
  ################################################################################
  #
  # survival plots
  #
  ################################################################################

  #
  # for (surv_os_dataset in ..) {
  #   output[['%s_plot']] = renderPlot({})
  # }
  surv_os_datasets <- covariate_mdata %>% 
    filter(phenotype=='overall_survival') %>% 
    pull(dataset)

  # for (acc in surv_os_datasets) {
  lapply(surv_os_datasets, function(acc) {
    output[[sprintf("%s_os_plot", acc)]] <- renderPlot({
      req(input$surv_os_tbl_rows_selected)

      df <- indiv_mdata[[acc]] %>%
        select(1, time=os_time, event=os_censor)

      sample_ids <- df %>%
        pull(1)

      gene_name <- surv_os_gene_selected()
      gene_expr <- as.numeric(expr_dat[gene_name, sample_ids])

      df$feature <- gene_expr 

      if (sum(!is.na(gene_expr)) > 0) {
        plot_survival(df, acc, "Overall Survival", gene_name, "Days", surv_expr_cutoffs, cfg$colors)
      }
    })
  })

  output$surv_os_plots <- renderUI({
    log_info("output$surv_os_plots")

    plts <- list()

    for (acc in surv_os_datasets) {
      plot_id <- sprintf("%s_os_plot", acc)
      plts <- c(plts, list(plotOutput(plot_id, height="360px")))
    }
    plts
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

      # tabPanel(
      #   "Disease Stage",
      #   fluidRow(
      #     column(
      #         width=6,
      #         withSpinner(DTOutput("disease_stage_tbl"))
      #     ),
      #     column(
      #         width=6,
      #         selectizeInput("gene_stg", "Gene:", choices=stage_rankings$Gene),
      #         helpText("Gene to visualize disease stage data for."),
      #         hr(),
      #         withSpinner(plotOutput("stage_plot", height="740px"))
      #     )
      #   )
      # ),
      tabPanel(
        "OS",
        fluidRow(
          column(
              width=6,
              withSpinner(DTOutput("surv_os_tbl"))
          ),
          column(
              width=6,
              tableOutput("surv_os_details_tbl"),
              withSpinner(htmlOutput("surv_os_plots"))
          ),
        ),
      ),
      # tabPanel(
      #   "Treatment",
      #   fluidRow(
      #     column(
      #         width=6,
      #         withSpinner(DTOutput("treatment_response_tbl"))
      #     ),
      #     column(
      #         width=6,
      #         fluidRow(
      #           column(
      #             width=4,
      #             selectizeInput("gene_trmt", "Gene:",
      #                            choices=treatment_response_rankings$Gene),
      #             helpText("Gene to visualize treatment response data for.")
      #           ),
      #           column(
      #            width=2,
      #            switchInput("show_indiv_treatment_plots", value=TRUE),
      #            helpText("Show individual plots.")
      #           )
      #         ),
      #         hr(),
      #         withSpinner(htmlOutput("treatment_response_plots"))
      #     )
      #   )
      # ),
      # tabPanel(
      #   "Cell Lines",
      #   fluidRow(
      #       column(
      #         width=4,
      #         selectizeInput("gene_hmcl", "Gene:",
      #                         choices=hmcl_expr$symbol),
      #         helpText("Gene to show cell line gene expression data for.")
      #       ),
      #       column(
      #         width=6,
      #         withSpinner(DTOutput("hmcl_expr_tbl"))
      #       ),
      #   )
      # ),
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
