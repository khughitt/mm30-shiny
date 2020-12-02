#
# MM25 Shiny UI
# KH Nov 2020
#
# TODO:
#
# - refactor..
# - simplify by renaming symbol/pathway to "feature_id"?
# - color gene symbols by pubmed counts, etc.?
#
library(annotables)
library(arrow)
library(DT)
library(ggdark)
library(ggforce)
library(ggrepel)
library(gridExtra)
library(heatmaply)
library(Hmisc)
library(plotly)
library(shinythemes)
library(shinycssloaders)
library(survival)
library(survminer)
library(tidyverse)
library(uwot)
library(yaml)

source("R/plotting.R")

options(stringsAsFactors = FALSE)
options(spinner.color="#00bc8c")
options(scipen = 2, digits = 3)
set.seed(1)

# ggplot theme
theme_dark <- function(base_font_size = 18) {
  dark_theme_gray(base_size = base_font_size) +
    theme(axis.text.x = element_text(angle = 90),
          legend.background = element_rect(fill = NA),
          plot.background = element_rect(fill = "#222222"),
          panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
          panel.grid.major = element_line(color = "#555555", size = 0.2),
          panel.grid.minor = element_line(color = "#555555", size = 0.2))
}

# heatmap theme
heatmap_theme <- dark_theme_gray() +
  theme(plot.background = element_rect(fill = "#222222"),
        panel.background = element_rect(fill = "#222222"),
        panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
        legend.background = element_rect(fill = "#222222"),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(color = "white"))

# feature levels
feature_levels <- c(Genes = "genes", Pathways = "gene_sets")

# plot colors
color_pal <- c("#36c764", "#c73699", "#3650c7", "#c7ac36", "#c7365f", "#bb36c7")

# options for survival plot upper/lower expression cutoffs
surv_expr_cutoffs <- seq(5, 50, by = 5)
names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")

#############################
#
# Shiny Server
#
#############################
server <- function(input, output, session) {
  #
  # Reactives
  #

  # load config
  cfg <- reactive({
    req(input$select_version)
    read_yaml(sprintf("config/config-%s.yml", input$select_version))
  })

  # load dataset and covariate metadata
  phenotype_metadata <- reactive({
    message("[phenotype_metadata]")

    dataset_metadata <- read_tsv(cfg()$dataset_metadata, col_types = cols())

    infile <- cfg()$phenotype_metadata
    message(sprintf("[phenotype_metadata] loading %s", infile))

    dat <- read_feather(infile)

    if ("feature_level" %in% colnames(dat)) {
      dat <- dat %>%
        filter(feature_level == "genes") %>%
        select(-feature_level)
    }

    # work-around: manually add MMRF to dataset metadata
    dataset_metadata <- rbind(dataset_metadata,
      c("MMRF", "MMRF CoMMpass Study IA14", NA, NA, NA, NA, NA, "GPL11154", NA,
        "https://research.themmrf.org/", "", ""))

    dataset_metadata <- dataset_metadata %>%
      rename(dataset = geo_id)

    dat <- dat %>%
      inner_join(dataset_metadata, by = "dataset")

    # add covariate id
    # dat$covariate_id <- sprintf("%s_%s_pval", dat$dataset, dat$phenotype)
    dat$covariate_id <- sprintf("%s_%s", dat$dataset, dat$phenotype)

    # add covariate category
    dat$category <- factor(dat$category)

    # add covariate cluster
    infile <- cfg()$phenotype_clusters
    message(sprintf("[phenotype_metadata] loading %s", infile))

    covariate_clusters <- read_feather(infile) %>%
      rename(covariate_id = covariate)

    dat <- dat %>%
      inner_join(covariate_clusters, by = 'covariate_id')

    # offset by one so that first cluster is "cluster 1" and convert to a factor
    dat$cluster <- factor(dat$cluster + 1)

    # get the number of genes and samples in processed expression data
    num_genes   <- lapply(gene_data(), nrow)
    num_samples <- lapply(gene_data(), ncol)

    dat$num_genes   <- as.numeric(num_genes[dat$dataset])
    dat$num_samples <- as.numeric(num_samples[dat$dataset])

    # add number of significant gene associations and fgsea results
    num_sig_assoc <- apply(mm25_feature_padjs(), 2, function (x) {
      sum(x < 0.01, na.rm = TRUE)
    })

    dat$num_sig_assoc <- num_sig_assoc[match(dat$covariate_id, names(num_sig_assoc))]

    # add number of significant fgsea terms
    fgsea_summary <- read_tsv(cfg()$fgsea_summary_indiv, col_types = cols())

    # backwards-compat (MM25 v1.0 - v1.1)
    fgsea_summary$field <- sub("_pval$", "", fgsea_summary$field)

    dat$num_sig_gsea <- fgsea_summary$num_sig[match(dat$covariate_id, fgsea_summary$field)]

    dat
  })

  # load MM25 combined pvals
  mm25_gene_pvals_combined <- eventReactive(input$select_version_subset, {
    infile <- cfg()$mm25_scores$genes[[input$select_version_subset]]

    message(sprintf("[mm25_gene_pvals_combined] loading %s", infile))

    read_feather(infile) %>%
      select(symbol, sumz_wt_pval, sumz_pval, sumlog_pval, min_pval, num_present, num_missing)
  })

  mm25_pathway_pvals_combined <- reactive({
    req(input$select_version_subset)

    infile <- cfg()$mm25_scores$pathways[[input$select_version_subset]]

    message(sprintf("[mm25_pathway_pvals_combined] loading %s", infile))

    dat <- read_feather(infile)

    # add links to msigdb pathway info pages
    msigdb_ids <- sub("[^_]*_", "", dat$gene_set)

    msigdb_links <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>",
                            msigdb_ids, dat$gene_set)
    dat %>%
      add_column(pathway = msigdb_links, .before = 1) %>%
      select(pathway, gene_set, sumz_wt_pval, sumz_pval, sumlog_pval, min_pval, num_present, num_missing)
  })

  mm25_gene_pvals_indiv <- reactive({
    infile <- cfg()$gene_pvals_indiv
    message(sprintf("[mm25_gene_pvals_indiv] loading %s", infile))
    read_feather(infile)
  })

  mm25_pathway_pvals_indiv <- reactive({
    infile <- cfg()$pathway_pvals_indiv
    message(sprintf("[mm25_pathway_pvals_indiv] loading %s", infile))
    read_feather(infile)
  })

  # individual dataset p-values
  mm25_feature_pvals <- reactive({
    message("[mm25_feature_pvals]")

    if (rv$plot_feature_level == 'genes') {
      dat <- mm25_gene_pvals_indiv() %>%
        rename(feature_id = symbol)
    } else {
      dat <- mm25_pathway_pvals_indiv() %>%
        rename(feature_id = gene_set)
    }

    # backwards-compat (v1.0 - v1.1)
    colnames(dat) <- sub("_pval$", "", colnames(dat))

    dat
  })

  mm25_feature_padjs <- reactive({
    message("[mm25_feature_padjs]")

    mm25_feature_pvals() %>%
      mutate_at(vars(-feature_id), p.adjust, method = "BH")
  })

  # list of covariates
  covariate_ids <- reactive({
    colnames(mm25_feature_pvals())[-1]
  })

  # load individual dataset configs;
  # TODO: store all needed paths, etc. in phenotype metadata file instead,
  # if possible
  dataset_cfgs <- reactive({
    cfg_infiles <- Sys.glob(file.path(cfg()$dataset_cfg_dir, "*.yml"))
    cfgs <- lapply(cfg_infiles, read_yaml)
    names(cfgs) <- sapply(cfgs, "[[", "name")

    cfgs
  })

  # load individual dataset gene expression data;
  # used to generate gene-level feature vs. phenotype plots
  gene_data <- reactive({
    message("[gene_data] init")

    gene_infiles <- lapply(dataset_cfgs(), function(x) {
      x$features$genes$rna
    })

    # for GEO datasets, make sure to use non-redundant (nr) versions of expression data
    mask <- startsWith(names(gene_infiles), "GSE")
    gene_infiles[mask] <- sub(".feather", "_nr.feather", gene_infiles[mask])

    # dat <- lapply(gene_infiles, read_feather)

    dat <- list()

    for (infile in gene_infiles) {
      message(sprintf("[gene_data] loading %s", infile))
      dat[[infile]] <- read_feather(infile)
    }
    names(dat) <- names(dataset_cfgs())

    dat
  })

  # list of all gene symbols
  all_genes <- reactive({
    sort(unique(unlist(lapply(gene_data(), '[', 'symbol'))))
  })

  # load individual dataset gene expression data;
  # used to generate gene-level feature vs. phenotype plots
  pathway_data <- reactive({
    message("[pathway_data] init")

    read_data <- function(x) {
      message(sprintf("[pathway_data] loading %s", x))

      if (endsWith(x, 'parquet')) {
        read_parquet(x)
      } else {
        read_feather(x)
      }
    }

    infiles <- lapply(dataset_cfgs(), function(x) {
      x$features$gene_sets$rna
    })

    dat <- lapply(infiles, read_data)
    names(dat) <- names(dataset_cfgs())

    dat
  })

  pheno_data <- reactive({
    dat <- list()

    for (dataset_id in names(dataset_cfgs())) {
      dataset_cfg <- dataset_cfgs()[[dataset_id]]
      infile <- dataset_cfg$phenotypes$path

      message(sprintf("[pheno_data] loading %s", infile))

      if (endsWith(infile, "tsv") || endsWith(infile, "tsv.gz")) {
        dat[[dataset_id]] <- read_tsv(infile, col_types = cols())
      } else if (endsWith(infile, "feather")) {
        dat[[dataset_id]] <- read_feather(infile)
      }
    }

    dat
  })

  selected_feature_pvals <- reactive({
    req(rv$plot_feature)

    mask <- mm25_feature_pvals()$feature_id == rv$plot_feature

    dat <- tibble(
      phenotype = covariate_ids(),
      pval = as.numeric(mm25_feature_pvals()[mask, -1]),
      padj = as.numeric(mm25_feature_padjs()[mask, -1])
    )

    dat <- dat %>%
      arrange(pval)
  })

  selected_phenotype_metadata <- reactive({
    req(rv$plot_covariate)

    # determine selected dataset / covariate
    dataset_id <- unlist(str_split(rv$plot_covariate, "_"))[1]
    pheno <- sub("_pval", "", sub(paste0(dataset_id, "_"), "", rv$plot_covariate))

    # retrieve relevant metadata
    phenotype_metadata() %>%
      filter(dataset ==  dataset_id & phenotype == pheno)
  })

  #
  # reactives - fgsea
  #
  fgsea_summary_indiv <- reactive({
    # dat <- read_tsv(cfg()$fgsea_summary_indiv, col_types = cols()) %>%
    #   arrange(desc(num_sig))
    message("[fgsea_summary_indiv]")

    phenotype_metadata() %>%
      select(dataset, phenotype, num_samples, num_sig_gsea, num_sig_assoc) %>%
      arrange(desc(num_sig_gsea))
  })

  fgsea_summary_combined <- reactive({
    req(input$select_fgsea_subset)

    infile <- cfg()$mm25_fgsea$summary[[input$select_fgsea_subset]]

    read_tsv(infile, col_types = cols()) %>%
      arrange(desc(num_sig))
  })

  fgsea_results_indiv <- reactive({
    message("[fgsea_results_indiv]")

    infile <- cfg()$fgsea_results_indiv

    message(sprintf("[fgsea_results_indiv] loading %s", infile))

    if (endsWith(infile, 'parquet')) {
      dat <- read_parquet(infile)
    } else {
      dat <- read_feather(infile)
    }

    # backwards-compat (v1.0 - v1.1)
    dat$field <- sub("_pval$", "", dat$field)

    dat %>%
        select(-dataset) %>%
        filter(padj < 0.05) %>%
        arrange(padj)
  })

  # fgsea_results_combined <- reactive({
  #   infile <- cfg()$fgsea_results_combined
  #
  #   if (endsWith(infile, 'parquet')) {
  #     dat <- read_parquet(infile)
  #   } else {
  #     dat <- read_feather(infile)
  #   }
  #
  #   dat %>%
  #     select(-dataset) %>%
  #     filter(padj < 0.05) %>%
  #     arrange(padj)
  # })

  fgsea_results_indiv_filtered <- reactive({
    req(input$select_fgsea_covariate)

    dat <- fgsea_results_indiv() %>%
        filter(field == input$select_fgsea_covariate) %>%
        select(-field)

    dat$pathway <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>",
                            dat$pathway, dat$pathway)

    dat
  })

  # fgsea_results_combined_filtered <- reactive({
  #   req(input$select_fgsea_ranking)
  #
  #   dat <- fgsea_results_combined() %>%
  #       filter(field == input$select_fgsea_ranking) %>%
  #       select(-field)
  #
  #   dat$pathway <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>",
  #                           dat$pathway, dat$pathway)
  #
  #   dat
  # })

  #
  # Text output
  #
  output$page_title <- renderText({
    req(input$select_version)
    sprintf("MM25 (%s)", input$select_version)
  })

  #
  # html output
  #
  output$covariate_summary <- renderUI({
    pheno_mdata <- selected_phenotype_metadata()

    tag_list <- tagList(
      tags$b(pheno_mdata$title),
      br(),
      tags$b(tags$a(href=pheno_mdata$url, target='_blank', pheno_mdata$dataset)),
      br()
    )

    if (!is.na(pheno_mdata$name)) {
      authors <- str_to_title(str_match(pheno_mdata$name, "[[:alpha:]]+"))
      year <- str_split(pheno_mdata$submission_date, " ", simplify = TRUE)[, 3]

      attribution <- sprintf("%s (%s)", authors, year)
      tag_list <- tagAppendChildren(tag_list, attribution, br())
    }

    tag_list <- tagAppendChildren(
      tag_list,
      "# Samples:",
      tags$b(pheno_mdata$num_samples),
      br(),
      "# Genes:",
      tags$b(pheno_mdata$num_genes),
      br()
    )

    if (!is.na(pheno_mdata$overall_design)) {
      tag_list <- tagAppendChildren(tag_list, hr(),
                                    tags$b("Overall Design:"), br(),
                                    pheno_mdata$abstract, br())
    }

    tags$div(tag_list)
  })

  #
  # Form fields
  #
  output$select_plot_covariate <- renderUI({
    req(rv$plot_feature)

    dat <- selected_feature_pvals()

    # select choices (e.g. "MMRF IA14 overall survival (p = 0.003)")
    opts <- selected_feature_pvals()$phenotype

    # remove dataset id prefix from covariate field
    # opts <- str_replace(opts, "[a-zA-Z0-9]+_", "")

    names(opts) <- sprintf("%s (padj = %0.3f)",
                           gsub("_", " ", sub("_pval", "", dat$phenotype)),
                           dat$padj)

    selectInput("select_plot_covariate", "Covariate:", choices = opts)
  })

  # output$select_version_subset <- renderUI({
  #   message("output$select_version_subset")
  #   opts <- names(cfg()$mm25_scores$genes)
  #   names(opts) <- Hmisc::capitalize(gsub('_', ' ', opts))
  #   selectInput("select_version_subset", "Subset:", choices = opts, selected = "all")
  # })

  output$fgsea_select_covariate <- renderUI({
    message("[fgsea_select_covariate]")

    fgsea_summary <- fgsea_summary_indiv() %>%
      select(dataset, phenotype, num_sig = num_sig_gsea) %>%
      mutate(covariate_id = sprintf("%s_%s", dataset, phenotype))
      # mutate(covariate_id = sprintf("%s_%s_pval", dataset, phenotype))

    opts <- covariate_ids()

    num_sig <- fgsea_summary$num_sig[match(opts, fgsea_summary$covariate_id)]

    # include number of significant terms in labels
    names(opts) <- sprintf("%s (#sig: %d)", opts, num_sig)

    # order labels in decreasing order of functional enrichment
    opts <- opts[match(fgsea_summary$covariate_id, opts)]

    selectInput("select_fgsea_covariate", "Covariate:", choices = opts)
  })

  # output$fgsea_select_ranking <- renderUI({
  #   fgsea_summary <- fgsea_summary_combined()
  #
  #   opts <- fgsea_summary$field
  #
  #   # include number of significant terms in labels
  #   names(opts) <- sprintf("%s (#sig: %d)", opts, fgsea_summary$num_sig)
  #
  #   # order labels in decreasing order of functional enrichment
  #   # opts <- opts[match(fgsea_summary$field, opts)]
  #
  #   selectInput("select_fgsea_ranking", "Ranking:", choices = opts)
  # })

  output$select_fgsea_subset <- renderUI({
    opts <- names(cfg()$mm25_fgsea$summary)
    names(opts) <- Hmisc::capitalize(gsub('_', ' ', opts))
    selectInput("select_fgsea_subset", "Subset (Affects Ranking Methods Only): ", choices = opts, selected = 'all')
  })

  # output$select_nlp_concept <- renderUI({
  #   selectInput("select_nlp_concept", "Disease Term", choices =
  #               colnames(pubtator_gene_disease())[-1], selected = "Multiple Myeloma")
  # })

  rv <- reactiveValues(
    # version_subset = NULL,
    plot_feature_level = NULL,
    plot_feature = NULL,
    plot_covariate = NULL
  )

  #
  # Manually setup / trigger select_version event
  #
  updateSelectInput(session, "select_version", "Version:",
                    choices = c("v2.0", "v2.1", "v3.0", "v3.1", "v3.2"), selected = "v3.2")

  # Using static set of GRCh38 genes from annotables; the datasets in MM25 also
  # include some other gene symbols, but there are generally quite rare and not
  # likely to be of interest
  updateSelectizeInput(session, "select_coex_gene1", choices = grch38$symbol,
                       selected = "MCL1", server = TRUE)

  updateSelectizeInput(session, "select_coex_gene2", choices = grch38$symbol,
                       selected = "PBXIP1", server = TRUE)

  #
  # Event Hanlders
  #

  # 1. version subset
  observeEvent(input$select_version, {
    req(cfg())

    # reset form fields when version changes
    message("observeEvent(input$select_version)")

    opts <- names(cfg()$mm25_scores$genes)
    names(opts) <- Hmisc::capitalize(gsub('_', ' ', opts))

    # in order to trigger the select version subset, first set it to NULL, then
    # repopulate the choices

    # if this doesn't work, try creating a (global) variable that indicates whether
    # the app has been reset, and if so, update each piece until ready..
    updateSelectInput(session, "select_version_subset", "Subset:",
                      choices = c("loading..."), selected = "loading...")

    updateSelectInput(session, "select_version_subset", "Subset:",
                      choices = opts, selected = "all")
  })

  observeEvent(input$select_version_subset, {
    req(input$select_version)

    if (input$select_version_subset != 'loading...') {
      rv$plot_feature_level <- "genes"

      message("observeEvent(input$select_version_subset)")

      updateSelectInput(session, "select_plot_feature_level", "Feature Level:",
                        choices = feature_levels, selected = "genes")
    }
  })

  observeEvent(input$select_plot_feature_level, {
    req(input$select_version_subset)

    rv$plot_feature <- NULL
    rv$plot_covariate <- NULL

    if (input$select_plot_feature_level != "") {
      rv$plot_feature_level <- input$select_plot_feature_level
    }

    message("observeEvent(input$select_plot_feature_level)")


    if (rv$plot_feature_level == "genes") {
      select_choices <- as.character(mm25_gene_pvals_combined()$symbol)
    } else {
      select_choices <- as.character(mm25_pathway_pvals_combined()$gene_set)
    }

    updateSelectizeInput(session, "select_plot_feature", choices = select_choices, server = TRUE)
  })


  observeEvent(input$select_plot_covariate, {
    message("observeEvent(input$select_plot_covariate)")
    rv$plot_covariate <- input$select_plot_covariate
  })

  observeEvent(input$select_plot_feature, {
    req(input$select_plot_feature_level)

    rv$plot_covariate <- NULL
    rv$plot_feature <- input$select_plot_feature

    message("observeEvent(input$select_plot_feature)")

    # dat <- selected_feature_pvals()
    #
    # # select choices (e.g. "MMRF IA14 overall survival (p = 0.003)")
    # opts <- selected_feature_pvals()$phenotype
    #
    # names(opts) <- sprintf("%s (padj = %0.3f)",
    #                        gsub("_", " ", sub("_pval", "", dat$phenotype)),
    #                        dat$padj)
    #
    # updateSelectInput(session, "select_plot_covariate", "Covariate:", choices = opts)
  })

  #
  # Tables
  #
  output$mm25_gene_pvals_combined_table <- renderDataTable({
    req(input$select_version_subset)
    req(input$select_table_format)

    float_cols <- paste0(c("min", "sumlog", "sumz", "sumz_wt"), "_pval")

    dat <- mm25_gene_pvals_combined()

    # convert to ranks, if requested
    if (input$select_table_format == "Ranks") {
      dat <- dat %>%
        mutate_at(vars(ends_with("_pval")), dense_rank)
    }

    # add biotype
    gene_biotypes <- grch38 %>%
      select(symbol, biotype) %>%
      group_by(symbol) %>%
      slice(1) %>%
      ungroup()

    dat <- dat %>%
      left_join(gene_biotypes, by = 'symbol') %>%
      select(symbol, biotype, everything())

    # construct data table
    out <- DT::datatable(dat, style = "bootstrap", options = list(pageLength = 15))

    if (input$select_table_format == "P-values") {
      out <- out %>%
        formatSignif(columns = float_cols, digits = 3)
    }

    out
  })

  output$mm25_pathway_pvals_combined_table <- renderDataTable({
    req(input$select_version_subset)
    # req(input$select_pathway_table_format)

    message("mm25_pathway_pvals_combined_table")

    float_cols <- paste0(c("min", "sumlog", "sumz", "sumz_wt"), "_pval")

    dat <- mm25_pathway_pvals_combined()

    # convert to ranks, if requested
    # if (input$select_pathway_table_format == "Ranks") {
    #   dat <- dat %>%
    #     mutate_at(vars(ends_with("_pval")), dense_rank)
    # }

    # construct data table
    DT::datatable(dat %>% select(-gene_set),
                         style = "bootstrap", escape = FALSE,
                         options = list(pageLength = 15)) %>%
      formatSignif(columns = float_cols, digits = 3)
  })

  output$fgsea_results_indiv <- renderDataTable({
    DT::datatable(fgsea_results_indiv_filtered(),
                  style = "bootstrap", escape = FALSE, options = list(pageLength = 15)) %>%
        formatSignif(columns = c("pval", "padj", "ES", "NES"), digits = 3)
  })

  # output$fgsea_results_combined <- renderDataTable({
  #   DT::datatable(fgsea_results_combined_filtered(),
  #                 style = "bootstrap", escape = FALSE, options = list(pageLength = 15)) %>%
  #       formatSignif(columns = c("pval", "padj", "ES", "NES"), digits = 5)
  # })

  output$fgsea_summary_table <- renderDataTable({
    req(input$select_fgsea_summary)

    if (input$select_fgsea_summary == "Covariates") {
        dat <- fgsea_summary_indiv()
    } else {
        req(input$select_fgsea_subset)
        dat <- fgsea_summary_combined()
    }

    DT::datatable(dat, style = "bootstrap", options = list(pageLength = 15))
  })

  output$nlp_pubtator_rankings <- renderDataTable({
    infile <- cfg()$pubtator$gene_disease

    message(sprintf("[nlp_pubtator_rankings] loading %s", infile))

    dat <- read_feather(infile) %>%
      rename(symbol = gene) %>%
      group_by(symbol) %>%
      summarize_all(max)

    gene_ranks <- mm25_gene_pvals_combined() %>%
      select(symbol, sumz_wt_pval) %>%
      mutate(rank = dense_rank(sumz_wt_pval)) %>%
      select(-sumz_wt_pval)

    dat <- gene_ranks %>%
      inner_join(dat, by = 'symbol')

    DT::datatable(dat, style = "bootstrap", options = list(pageLength = 15))
  })

  #
  # Plots
  #
  output$feature_plot <- renderPlot({
    req(rv$plot_covariate)
    req(input$select_survival_expr_cutoffs)

    message("feature_plot")

    pheno <- rv$plot_covariate

    # backwards-compatibility
    pheno <- sub("_pval", "", pheno)

    # get dataset and covariate names
    dataset_id <- unlist(str_split(pheno, "_"))[1]
    covariate <- sub(paste0(dataset_id, "_"), "", pheno)

    # feature data
    if (rv$plot_feature_level == 'genes') {
      feat_dat <- gene_data()[[dataset_id]] %>%
        filter(symbol == rv$plot_feature) %>%
        select(-symbol) %>%
        as.numeric()
    } else {
      feat_dat <- pathway_data()[[dataset_id]] %>%
        filter(gene_set == rv$plot_feature) %>%
        select(-gene_set) %>%
        as.numeric()
    }

    # get config for selected feature-phenotype association
    assoc_cfg <- dataset_cfgs()[[dataset_id]]$phenotypes$associations[[covariate]]
    assoc_method <- assoc_cfg$method

    # assoc_method <- phenotype_metadata() %>%
    #   filter(dataset == dataset_id & phenotype == covariate) %>%
    #   pull(method)

    if (assoc_method == "survival") {
      # survival plot
      pheno_dat <- pheno_data()[[dataset_id]]

      time_dat <- pull(pheno_dat, assoc_cfg$params$time)
      event_dat <- pull(pheno_dat, assoc_cfg$params$event)

      dat <- data.frame(feature = feat_dat, time = time_dat, event = event_dat)

      plot_survival(dat, dataset_id, covariate, rv$plot_feature,
                    input$select_survival_expr_cutoffs, color_pal, theme_dark)
    } else {
      # violin plot
      response <- factor(pull(pheno_data()[[dataset_id]], assoc_cfg$params$field))

      dat <- data.frame(feature = feat_dat, response)

      plot_categorical(dat, dataset_id, covariate, rv$plot_feature, color_pal, theme_dark)
    }
  })

  output$cov_similarity_plot <- renderPlotly({
    req(input$select_cov_similarity_feat_type)
    req(input$select_cov_similarity_cor_method)

    # determine dataset to use
    if (input$select_cov_similarity_feat_type == "Genes") {
      dat <- mm25_gene_pvals_indiv()
    } else {
      dat <- mm25_pathway_pvals_indiv()

    }

    # drop gene names; not needed here
    dat <- dat[, -1]

    # determine covariate functional groups from labels;
    # TODO: load dataset metadata table including groups and use that instead
    cnames <- sub("_pval", "", colnames(dat))

    dataset_ids <- paste(phenotype_metadata()$dataset,
                         phenotype_metadata()$phenotype, sep="_")
    cov_categories <- phenotype_metadata()$category[match(cnames, dataset_ids)]
    cov_clusters   <- phenotype_metadata()$cluster[match(cnames, dataset_ids)]
    annot_row <- data.frame(type = factor(cov_categories))
    annot_col <- data.frame(cluster = cov_clusters)

    # cov_categories <- rep("", ncol(dat) - 1)
    # cov_categories[grepl("survival|died", cnames)] <- "survival"
    # cov_categories[grepl("treatment|response", cnames)] <- "treatment"
    # cov_categories[grepl("status|stage|ecog|pfs_event|relapsed", cnames)] <- "stage"

    # cov_methods <- phenotype_metadata()$method[match(cnames, dataset_ids)]
    # annot_col <- data.frame(method = factor(cov_methods))

    # -log10 transform p-values (clipping at 1-E20)
    dat[dat < 1E-20] <- 1E-20
    dat <- -log10(dat)

    # dataset category color palette
    category_pal <- color_pal[1:3]
    names(category_pal) <- c("survival", "treatment", "disease_stage")

    # cluster shapes to use
    cluster_shapes <- c(1, 15, 2, 16, 17)[1:length(unique(cov_clusters))]

    #
    # generate similarity plot
    #

    # 1. heatmap
    if (input$select_cov_similarity_plot_type == 'Heatmap') {
      # generate covariate correlation matrix
      # TODO; make -log10 transform optional..
      cor_method <- tolower(input$select_cov_similarity_cor_method)
      cor_mat <- cor(dat, method = cor_method, use = "pairwise.complete.obs")


      # render heatmap
      heatmaply(cor_mat, row_side_colors = annot_row, row_side_palette = category_pal,
                col_side_colors = annot_col, width = 1200, height = 1400,
                heatmap_layers = heatmap_theme,
                dendrogram_layers = list(
                  scale_color_manual(values = c("#b2b2b2", "#b2b2b2")),
                  heatmap_theme,
                  theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
                ),
                side_color_layers = heatmap_theme)
    } else if (input$select_cov_similarity_plot_type == 'PCA') {
      # 2. PCA plot

      # first, drop rows with missing values
      dat <- dat[complete.cases(dat), ]

      pca <- prcomp(t(dat), scale = TRUE)$x[, 1:2]
      colnames(pca) <- c("PC1", "PC2")

      pca <- cbind(annot_row, annot_col, pca)
      pca$covariate <- rownames(pca)

      ggplot(pca, aes(x = PC1, y = PC2, color = type, shape = cluster, label = covariate)) +
        geom_point(size = 3) +
        geom_text_repel() +
        ggtitle("MM25 Covariate PCA Plot") +
        scale_shape_manual(values = cluster_shapes) +
        scale_color_manual(values = category_pal) +
        theme_dark()
    } else if (input$select_cov_similarity_plot_type == 'UMAP') {
      # 3. UMAP plot
      set.seed(1)

      # first, drop rows with missing values
      dat <- dat[complete.cases(dat), ]

      cov_umap <- umap(t(dat), n_neighbors = 8, scale = FALSE)

      colnames(cov_umap) <- c("UMAP1", "UMAP2")

      cov_umap <- cbind(annot_row, annot_col, cov_umap)
      cov_umap$covariate <- colnames(dat)

      ggplot(cov_umap, aes(x = UMAP1, y = UMAP2, color = type, shape = cluster, label = covariate)) +
        geom_point(size = 3) +
        geom_text_repel() +
        ggtitle("MM25 Covariate UMAP Plot") +
        scale_shape_manual(values = cluster_shapes) +
        scale_color_manual(values = category_pal) +
        theme_dark()
    }
  })

  output$feature_pval_dists <- renderPlot({
    req(input$select_feature_pval_dist_type)

    set.seed(1)

    # adjusted / unadjusted p-values
    if (input$select_feature_pval_dist_type == "Adjusted (BH)") {
        dat <- mm25_feature_padjs()  %>%
        sample_n(1000) %>%
        pivot_longer(-feature_id, names_to = "covariate", values_to ="pvalue")
    } else {
        dat <- mm25_feature_pvals() %>%
        sample_n(1000) %>%
        pivot_longer(-feature_id, names_to = "covariate", values_to ="pvalue")
    }

    dat$dataset <- str_split(dat$covariate, "_", simplify = TRUE)[, 1]
    # dat$dataset[dat$dataset == "MMRF"] <- "MMRF_IA14"
    dat$dataset <- factor(dat$dataset)

    dat$covariate <- factor(sub("_pval", "",  dat$covariate))

    dat$label <- trimws(gsub("_", " ", str_replace(dat$covariate, as.character(dat$dataset), "")))
    dat$label <- sprintf("%s (%s)", dat$label, dat$dataset)

    dat$method <- "Logit"
    dat$method[grepl("survival", dat$covariate)] <- "Survival"
    dat$method <- factor(dat$method)

    npages <- ceiling((ncol(mm25_feature_pvals()) - 1) / 8)

    plts <- list()

    for (i in 1:npages) {
      plts[[i]] <- ggplot(dat, aes(x = pvalue, fill = method)) +
        geom_density(alpha = 0.8) +
        facet_wrap_paginate(~label, ncol = 2, nrow = 4, scales = "free_y", page = i) +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        theme_dark()
    }

    print(grid.arrange(grobs = plts, ncol = 1))
  })

  output$gene_coex_plots <- renderPlot({
    req(gene_data())
    req(input$select_coex_gene1)
    req(input$select_coex_gene2)

    gene1 <- input$select_coex_gene1
    gene2 <- input$select_coex_gene2

    # iterate over datasets and generate a plot for each one including both genes
    plts <- list()

    for (dat_name in names(gene_data())) {
      dat <- gene_data()[[dat_name]] %>%
        filter(symbol %in% c(gene1, gene2))

      if (nrow(dat) == 2) {
        dat_long <- dat %>%
          column_to_rownames('symbol') %>%
          t() %>%
          as.data.frame()

        # compute pearson correlation
        coex_cor <- cor(dat_long[, 1], dat_long[, 2], use = 'pairwise.complete')

        plts[[dat_name]] <- ggplot(dat_long, aes_string(gene1, gene2)) +
            geom_point() +
            geom_smooth(method = lm) +
            ggtitle(sprintf("%s (cor = %0.2f)", dat_name, coex_cor)) +
            theme_dark()
      }
    }

    print(grid.arrange(grobs = plts, ncol = 2))
  })

  output$fgsea_summary_plot <- renderPlotly({
    message("[fgsea_summary_plot]")

    # retrieve indiv and combined fgsea results
    indiv_dat <- phenotype_metadata() %>%
      select(field = covariate_id, num_sig = num_sig_gsea)

    dat <- rbind(cbind(indiv_dat, type = "individual"),
                 cbind(fgsea_summary_combined(), type = "combined"))

    dat$type <- factor(dat$type)

    sorted_fields <- dat %>%
      arrange(num_sig) %>%
      pull(field)

    # display in order of increasing significance
    dat$field <- factor(dat$field, levels = sorted_fields)

    # scale_color_manual(values = color_pal) +
    ggplot(dat, aes(x = field, y = num_sig, fill = type)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = color_pal) +
      xlab("Gene Ranking") +
      ylab("# Significant GSEA terms") +
      theme_dark(14)
  })
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
      id = "main",
      theme = shinytheme("darkly"),
      title = textOutput("page_title"),
      windowTitle = "MM25",

      tabPanel(
        "Genes",
        withSpinner(dataTableOutput("mm25_gene_pvals_combined_table"))
      ),
      tabPanel(
        "Pathways",
        withSpinner(dataTableOutput("mm25_pathway_pvals_combined_table"))
      ),
      tabPanel(
        "Visualize",
        fluidRow(
          column(
            width = 4,
            selectInput("select_plot_feature_level", "Feature Level:", choices = NULL),
            helpText("Feature level to visualize."),
            hr(),
            selectizeInput("select_plot_feature", "Feature:", choices = NULL),
            helpText("Feature to visualize."),
            hr(),
            uiOutput("select_plot_covariate"),
            helpText("Phenotype / covariate to compare feature expression or SNP counts against."),
            hr(),
            selectInput("select_survival_expr_cutoffs", "Upper/Lower Feature Expression Cutoffs:",
                        choices=surv_expr_cutoffs, selected = 25),
            helpText("Upper and lower feature expression percentile cutoffs to use for two survival groups."),
            hr(),
            uiOutput("covariate_summary")
          ),
          column(
            width = 8,
            withSpinner(plotOutput("feature_plot", height = "760px"))
          )
        )
      ),
      navbarMenu(
        "Covariates",
        tabPanel(
          "P-value Distributions",
          selectInput("select_feature_pval_dist_type", "P-value type:",
                      choices = c("Adjusted (BH)", "Unadjusted"), selected = "Unadjusted"),
          withSpinner(plotOutput("feature_pval_dists", height = "4000px"))
        ),
        tabPanel(
          "Covariate Similarity Plot",
          column(
            width = 3,
            selectInput("select_cov_similarity_plot_type", "Plot type:",
                        choices = c("Heatmap", "PCA", "UMAP"), selected = "PCA"),
            selectInput("select_cov_similarity_feat_type", "Feature type:",
                        choices = c("Genes", "Pathways"), selected = "Genes"),
            selectInput("select_cov_similarity_cor_method", "Correlation type:",
                        choices = c("Pearson", "Spearman"), selected = "Pearson")
          ),
          column(
            width = 9,
            withSpinner(plotlyOutput("cov_similarity_plot", height = "800px"))
          )
        )
      ),
      tabPanel(
        "Co-expression",
        column(
          width = 3,
          selectizeInput("select_coex_gene1", "Gene 1", choices = NULL),
          selectizeInput("select_coex_gene2", "Gene 2", choices = NULL)
        ),
        column(
          width = 9,
          withSpinner(plotOutput("gene_coex_plots", height = "4000px"))
        ),
      ),
      tabPanel(
        # "NLP",
        # column(
        #   width = 3,
        #   uiOutput("select_nlp_concept")
        # )
        "NLP",
        column(
          width = 12,
          withSpinner(dataTableOutput("nlp_pubtator_rankings"))
        ),
      ),
      navbarMenu(
        "Functional Enrichment",
        tabPanel(
          "Phenotypes (P-values)",
          uiOutput("fgsea_select_covariate"),
          withSpinner(dataTableOutput("fgsea_results_indiv"))
        ),

        #
        # May 27, 2020:
        #
        # Disabling this section for now until refactoring can be performed and it can
        # be more easily extended to support MM25 category/cluster subsets.
        #
        # tabPanel(
        #   "MM25 Rankings",
        #   # uiOutput("fgsea_select_ranking"),
        #   fluidRow(
        #     column(width = 2, uiOutput("select_fgsea_subset")),
        #     column(width = 2, uiOutput("fgsea_select_ranking"))
        #   ),
        #   withSpinner(dataTableOutput("fgsea_results_combined"))
        # ),
        tabPanel(
          "Summary",
          fluidRow(
            column(width = 3,
                   selectInput("select_fgsea_summary", "Display: ",
                               choices = c("Covariates", "Ranking Methods"))),
            column(width = 3, uiOutput("select_fgsea_subset"))
          ),
          fluidRow(
            column(
              width = 6,
              withSpinner(dataTableOutput("fgsea_summary_table"))
            ),
            column(
              width = 6,
              withSpinner(plotlyOutput("fgsea_summary_plot", height = "800px"))
            )
          )
        )
      ),
      tabPanel(
        "Settings",
        selectInput("select_version", "Version:", choices = NULL),
        selectInput("select_version_subset", "Subset:", choices = NULL),
        selectInput("select_table_format", "Display:", choices = c("P-values", "Ranks"), selected = "Ranks")
      )
    )
  )
}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
