#
# Feature weights shiny UI
# KH Jan 2020
#
library(arrow)
library(DT)
library(feather)
library(ggdark)
library(ggforce)
library(gridExtra)
library(heatmaply)
library(plotly)
library(shinythemes)
library(shinycssloaders)
library(survival)
library(survminer)
library(tidyverse)
library(yaml)

options(stringsAsFactors = FALSE)
options(spinner.color="#00bc8c")
options(digits = 3)
set.seed(1)

# ggplot theme
theme_dark <- dark_theme_gray(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90),
        legend.background = element_rect(fill = NA),
        plot.background = element_rect(fill = "#222222"),
        panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
        panel.grid.major = element_line(color = "#555555", size = 0.2),
        panel.grid.minor = element_line(color = "#555555", size = 0.2))

# heatmap theme
heatmap_theme <- dark_theme_gray() +
  theme(plot.background = element_rect(fill = "#222222"),
        panel.background = element_rect(fill = "#222222"),
        panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
        legend.background = element_rect(fill = "#222222"),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(color = "white"))

# plot colors
color_pal <- c('#36c764', '#c73699', '#3650c7', '#c7ac36', '#c7365f', '#bb36c7')

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
  covariate_metadata <- reactive({
    dataset_metadata <- read_tsv(cfg()$dataset_metadata, col_types = cols())

    dat <- read_feather(cfg()$covariate_metadata)

    # work-around: manually add MMRF to dataset metadata
    dataset_metadata <- rbind(dataset_metadata,
      c('MMRF', 'MMRF CoMMpass Study IA14', NA, NA, NA, NA, NA, "GPL11154", NA,
        "https://research.themmrf.org/", "", ""))

    dataset_metadata <- dataset_metadata %>%
      rename(dataset = geo_id)

    dat <- dat %>%
      inner_join(dataset_metadata, by = 'dataset')

    dat$category <- factor(dat$category)

    # get the number of genes and samples in processed expression data
    num_genes   <- lapply(gene_data(), nrow)
    num_samples <- lapply(gene_data(), ncol)

    dat$num_genes   <- as.numeric(num_genes[dat$dataset])
    dat$num_samples <- as.numeric(num_samples[dat$dataset])

    dat
  })

  # load MM25 combined pvals

  # genes
  # mm25_gene_scores <- read_feather(cfg()$gene_scores)

  mm25_genes <- reactive({
    read_feather(cfg()$gene_pvals)
  })

  # pathways
  # mm25_pathway_scores <- read_feather(cfg()$pathway_scores)

  mm25_pathways <- reactive({
    dat <- read_feather(cfg()$pathway_pvals)

    # add links to msigdb pathway info pages
    msigdb_ids <- sub("[^_]*_", "", dat$gene_set)

    msigdb_links <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>",
                            msigdb_ids, dat$gene_set)
    dat <- dat %>%
      add_column(pathway = msigdb_links, .before = 1)

    dat
  })
    # inner_join(mm25_pathway_scores, by = 'gene_set')

  # clean-up
  # rm(mm25_gene_scores, mm25_pathway_scores)

  # individual dataset p-values
  mm25_gene_pvals <- reactive({
    read_feather(cfg()$gene_pvals_indiv)
  })

  mm25_pathway_pvals <- reactive({
    read_feather(cfg()$pathway_pvals_indiv)
  })

  mm25_gene_padjs <- reactive({
    mm25_gene_pvals() %>%
      mutate_at(vars(-symbol), p.adjust, method = 'BH')
  })

  # list of covariates
  covariate_ids <- reactive({
    colnames(mm25_gene_pvals())[-1]
  })

  # load individual dataset configs;
  # TODO: store all needed paths, etc. in phenotype metadata file instead,
  # if possible
  dataset_cfgs <- reactive({
    cfg_infiles <- Sys.glob(file.path(cfg()$dataset_cfg_dir, '*.yml'))
    cfgs <- lapply(cfg_infiles, read_yaml)
    names(cfgs) <- sapply(cfgs, '[[', 'name')

    cfgs
  })

  # load individual dataset gene expression data;
  # used to generate gene-level feature vs. phenotype plots
  gene_data <- reactive({
    gene_infiles <- lapply(dataset_cfgs(), function(x) {
      x$features$genes$rna
    })

    # for GEO datasets, make sure to use non-redundant (nr) versions of expression data
    mask <- startsWith(names(gene_infiles), "GSE")
    gene_infiles[mask] <- sub('.feather', '_nr.feather', gene_infiles[mask])

    dat <- lapply(gene_infiles, read_feather)
    names(dat) <- names(dataset_cfgs())

    dat
  })

  pheno_data <- reactive({
    dat <- list()

    for (dataset in names(dataset_cfgs())) {
      dataset_cfg <- dataset_cfgs()[[dataset]]

      infile <- dataset_cfg$phenotypes$path

      if (endsWith(infile, 'tsv') || endsWith(infile, 'tsv.gz')) {
        dat[[dataset]] <- read_tsv(infile, col_types = cols())
      } else if (endsWith(infile, 'feather')) {
        dat[[dataset]] <- read_feather(infile)
      }
    }

    dat
  })

################################################

  gene_pvals <- reactive({
    req(input$select_gene)

    mask <- mm25_gene_pvals()$symbol == input$select_gene

    dat <- tibble(
      phenotype = covariate_ids(),
      pval = as.numeric(mm25_gene_pvals()[mask, -1]),
      padj = as.numeric(mm25_gene_padjs()[mask, -1])
    )

    dat <- dat %>%
      arrange(pval)
  })

  selected_covariate_metadata <- reactive({
    req(input$select_covariate)

    # determine selected dataset / covariate
    dataset_id <- unlist(str_split(input$select_covariate, '_'))[1]
    covariate_id <- sub('_pval', '', sub(paste0(dataset_id, '_'), '', input$select_covariate))

    # retrieve relevant metadata
    covariate_metadata() %>%
      filter(dataset ==  dataset_id & phenotype == covariate_id)
  })

  # fgsea_summary_indiv <- reactive({
  #   res <- read_tsv(cfg()$fgsea_pvals_indiv, col_types = cols())
  #   res %>%
  #     arrange(desc(num_sig))
  # })
  #
  # fgsea_summary_combined <- reactive({
  #   pvals <- read_tsv(cfg()$fgsea_pvals_combined, col_types = cols())
  #   scores <- read_tsv(cfg()$fgsea_scores_combined, col_types = cols())
  #
  #   rbind(pvals, scores) %>%
  #     arrange(desc(num_sig))
  # })
  #
  # fgsea_results_indiv <- reactive({
  #   read_parquet(cfg()$fgsea_results_indiv) %>%
  #       select(-dataset) %>%
  #       filter(padj < 0.05) %>%
  #       arrange(padj)
  # })
  #
  # fgsea_results_combined <- reactive({
  #   pvals <- read_parquet(cfg()$fgsea_results_pvals) %>%
  #     select(-dataset) %>%
  #     filter(padj < 0.05) %>%
  #     arrange(padj)
  #
  #   scores <- read_parquet(cfg()$fgsea_results_scores) %>%
  #     select(-dataset) %>%
  #     filter(padj < 0.05) %>%
  #     arrange(padj)
  #
  #   rbind(pvals, scores)
  # })
  #
  # fgsea_results_indiv_filtered <- reactive({
  #   req(input$select_fgsea_covariate)
  #
  #   dat <- fgsea_results_indiv() %>%
  #       filter(field == input$select_fgsea_covariate) %>%
  #       select(-field)
  #
  #   dat$pathway <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>",
  #                           dat$pathway, dat$pathway)
  #
  #   dat
  # })
  #
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
  # HTML output
  #
  output$covariate_summary <- renderUI({
    pheno_mdata <- selected_covariate_metadata()

    tag_list <- tagList(
      tags$b(pheno_mdata$title),
      br(),
      tags$b(tags$a(href=pheno_mdata$url, target='_blank', pheno_mdata$dataset)),
      br()
    )

    if (!is.na(pheno_mdata$name)) {
      authors <- str_to_title(str_match(pheno_mdata$name, '[[:alpha:]]+'))
      year <- str_split(pheno_mdata$submission_date, ' ', simplify = TRUE)[, 3]

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
  output$select_covariate <- renderUI({
    req(input$select_gene)

    dat <- gene_pvals()

    # select choices (e.g. "MMRF IA14 overall survival (p = 0.003)")
    opts <- gene_pvals()$phenotype

    names(opts) <- sprintf("%s (padj = %0.3f)",
                           gsub("_", " ", sub("_pval", "", dat$phenotype)),
                           dat$padj)

    selectInput('select_covariate', "Covariate:", choices = opts)
  })

  # output$fgsea_select_covariate <- renderUI({
  #   fgsea_summary <- fgsea_summary_indiv()
  #
  #   opts <- covariate_ids()
  #
  #   num_sig <- fgsea_summary$num_sig[match(opts, fgsea_summary$field)]
  #
  #   # include number of significant terms in labels
  #   names(opts) <- sprintf("%s (#sig: %d)", opts, num_sig)
  #
  #   # order labels in decreasing order of functional enrichment
  #   opts <- opts[match(fgsea_summary$field, opts)]
  #
  #   selectInput("select_fgsea_covariate", "Covariate:", choices = opts)
  # })
  #
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

  observeEvent(input$select_version, {
    updateSelectizeInput(session, 'select_gene',
                         choices = as.character(mm25_genes()$symbol), server = TRUE)
  })

  #
  # Tables
  #
  output$mm25_genes <- renderDataTable({
    req(input$select_gene_table_format)

    float_cols <- paste0(c('mean', 'median', 'min', 'sumlog', 'sumz'), '_pval')

    dat <- mm25_genes()

    # convert to ranks, if requested
    if (input$select_gene_table_format == 'Ranks') {
      dat <- dat %>%
        mutate_at(vars(ends_with("_pval")), dense_rank)
    }

    # construct data table
    out <- DT::datatable(dat, style = 'bootstrap', options = list(pageLength = 15))

    if (input$select_gene_table_format == 'P-values') {
      out <- out %>% 
        formatRound(columns = float_cols, digits = 5)
    }

    out
  })
  output$mm25_pathways <- renderDataTable({
    req(input$select_pathway_table_format)

    float_cols <- paste0(c('mean', 'median', 'min', 'sumlog', 'sumz'), '_pval')

    dat <- mm25_pathways()

    # convert to ranks, if requested
    if (input$select_pathway_table_format == 'Ranks') {
      dat <- dat %>%
        mutate_at(vars(ends_with("_pval")), dense_rank)
    }

    # construct data table
    out <- DT::datatable(dat %>% select(-gene_set),
                  style = 'bootstrap', escape = FALSE, options = list(pageLength = 15))

    if (input$select_pathway_table_format == 'P-values') {
      out <- out %>%
        formatRound(columns = float_cols, digits = 5)
    }

    out
  })

  output$fgsea_results_indiv <- renderDataTable({
    DT::datatable(fgsea_results_indiv_filtered(),
                  style = 'bootstrap', escape = FALSE, options = list(pageLength = 15)) %>%
        formatRound(columns = c('pval', 'padj', 'ES', 'NES'), digits = 5)
  })

  output$fgsea_results_combined <- renderDataTable({
    DT::datatable(fgsea_results_combined_filtered(),
                  style = 'bootstrap', escape = FALSE, options = list(pageLength = 15)) %>%
        formatRound(columns = c('pval', 'padj', 'ES', 'NES'), digits = 5)
  })

  output$fgsea_summary_output <- renderDataTable({
    req(input$select_fgsea_summary)

    if (input$select_fgsea_summary == "Covariates") {
        dat <- fgsea_summary_indiv()
    } else {
        dat <- fgsea_summary_combined()
    }

    DT::datatable(dat, style = 'bootstrap', options = list(pageLength = 15))
  })

  #
  # Plots
  #
  output$gene_plot <- renderPlot({
    req(input$select_covariate)
    pheno <- input$select_covariate

    # work-around: detect MMRF prefix (includes underscore)
    # parts <- unlist(str_split(sub('_pval', '', pheno), '_'))
    pheno <- sub('_pval', '', pheno)

    dataset <- unlist(str_split(pheno, '_'))[1]
    covariate <- sub(paste0(dataset, '_'), '', pheno)

    # config
    assoc_cfg <- dataset_cfgs()[[dataset]]$phenotypes$associations[[covariate]]

    if (assoc_cfg$method == 'survival') {
      # survival data
      feature <- gene_data()[[dataset]] %>%
        filter(symbol == input$select_gene) %>%
        select(-symbol) %>%
        as.numeric()

      time_dat <- pull(pheno_data()[[dataset]], assoc_cfg$params$time)
      event_dat <- pull(pheno_data()[[dataset]], assoc_cfg$params$event)

      dat <- data.frame(feature, time = time_dat, event = event_dat)

      # divide gene expression into quartiles
      expr_quartiles <- quantile(dat$feature)

      dat$Expression <- ""
      dat$Expression[dat$feature <= expr_quartiles[2]] <- "Lower 25%"
      dat$Expression[dat$feature >= expr_quartiles[4]] <- "Upper 25%"

      # drop all data except for upper and lower quartiles
      dat <- dat %>%
        filter(Expression != "")

      dat$Expression <- factor(dat$Expression)

      cov_label <- str_to_title(gsub('_', ' ', covariate))
      plt_title <- sprintf("%s: %s vs. %s", dataset, cov_label, input$select_gene)

      # perform fit on binarized data
      fit <- survfit(Surv(time, event) ~ Expression, data = dat)

      # display a kaplan meier plot for result
      ggsurvplot(fit, data = dat, ggtheme = theme_dark, palette = color_pal,
                 title = plt_title,
                 legend = 'bottom', legend.title = 'Legend')
    } else {
      # logistic regression fit plot
      feature <- gene_data()[[dataset]] %>%
        filter(symbol == input$select_gene) %>%
        select(-symbol) %>%
        as.numeric()

      # pheno data
      response <- pull(pheno_data()[[dataset]], assoc_cfg$params$field)

      dat <- data.frame(feature, response)

      dat$response <- factor(dat$response)

      # drop any entries with missing values
      dat <- dat[!is.na(dat$response), ]

      # draw violin + jitter plot
      set.seed(1)

      ggplot(dat, aes(x = response, y = feature)) +
        # geom_boxplot(aes(fill = response, color = response), outlier.shape = NA) +
        geom_violin(aes(fill = response, color = response), alpha = 0.5, draw_quantiles = c(0.5)) +
        geom_jitter(aes(color = response), alpha = 0.8) +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        ggtitle(sprintf("%s: %s vs. %s", dataset, input$select_gene, covariate)) +
        xlab(covariate) +
        ylab(sprintf("%s expression", input$select_gene)) +
        theme_dark
        # dark_theme_gray(base_size = 18) +
        # theme(axis.text.x = element_text(angle = 90),
        #       plot.background = element_rect(fill = "#222222"),
        #       panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
        #       panel.grid.major = element_line(color = "#555555", size = 0.2),
        #       panel.grid.minor = element_line(color = "#555555", size = 0.2))
    }
  })

  output$cov_similarity_heatmap <- renderPlotly({
    req(input$select_cov_similarity_feat_type)
    req(input$select_cov_similarity_cor_method)

    # determine dataset to use
    if (input$select_cov_similarity_feat_type == 'Genes') {
      dat <- mm25_gene_pvals()
    } else {
      dat <- mm25_pathway_pvals()
    }

    # determine covariate functional groups from labels;
    # TODO: load dataset metadata table including groups and use that instead
    cnames <- sub('_pval', '', colnames(dat)[-1])

    dataset_ids <- paste(covariate_metadata()$dataset,
                         covariate_metadata()$phenotype, sep='_')
    cov_categories <- covariate_metadata()$category[match(cnames, dataset_ids)]
    annot_row <- data.frame(type = factor(cov_categories))

    # cov_categories <- rep("", ncol(dat) - 1)
    # cov_categories[grepl('survival|died', cnames)] <- 'survival'
    # cov_categories[grepl('treatment|response', cnames)] <- 'treatment'
    # cov_categories[grepl('status|stage|ecog|pfs_event|relapsed', cnames)] <- 'stage'

    # cov_methods <- covariate_metadata()$method[match(cnames, dataset_ids)]
    # annot_col <- data.frame(method = factor(cov_methods))

    # generate covariate correlation matrix
    # TODO; make -log10 transform optional..
    cor_method <- tolower(input$select_cov_similarity_cor_method)
    cor_mat <- cor(-log10(pmax(dat[, -1], 1E-20)), method = cor_method, use = 'pairwise.complete.obs')

    # row annotations
    row_pal <- color_pal[1:3]
    names(row_pal) <- c('survival', 'treatment', 'disease_stage')

    # render heatmap
    heatmaply(cor_mat, row_side_colors = annot_row, row_side_palette = row_pal,
              heatmap_layers = heatmap_theme,
              dendrogram_layers = list(
                scale_color_manual(values = c('#b2b2b2', '#b2b2b2')),
                heatmap_theme,
                theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
              ),
              side_color_layers = heatmap_theme)
  })

  output$gene_pval_dists <- renderPlot({
    req(input$select_gene_pval_dist_type)

    set.seed(1)

    # adjusted / unadjusted p-values
    if (input$select_gene_pval_dist_type == "Adjusted (BH)") {
        dat <- mm25_gene_padjs()  %>%
        sample_n(1000) %>%
        pivot_longer(-symbol, names_to = 'covariate', values_to ='pvalue')
    } else {
        dat <- mm25_gene_pvals() %>%
        sample_n(1000) %>%
        pivot_longer(-symbol, names_to = 'covariate', values_to ='pvalue')
    }

    dat$dataset <- str_split(dat$covariate, '_', simplify = TRUE)[, 1]
    # dat$dataset[dat$dataset == 'MMRF'] <- 'MMRF_IA14'
    dat$dataset <- factor(dat$dataset)

    dat$covariate <- factor(sub('_pval', '',  dat$covariate))

    dat$label <- trimws(gsub("_", " ", str_replace(dat$covariate, as.character(dat$dataset), '')))
    dat$label <- sprintf("%s (%s)", dat$label, dat$dataset)

    dat$method <- 'Logit'
    dat$method[grepl('survival', dat$covariate)] <- 'Survival'
    dat$method <- factor(dat$method)

    npages <- ceiling((ncol(mm25_gene_pvals()) - 1) / 9)

    plts <- list()

    for (i in 1:npages) {
      plts[[i]] <- ggplot(dat, aes(x = pvalue, fill = method)) +
        geom_density(alpha = 0.8) +
        facet_wrap_paginate(~label, ncol = 3, nrow = 3, scales = 'free_y', page = i) +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        theme_dark
    }

    print(grid.arrange(grobs = plts, ncol = 1))
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
      windowTitle = "MM25 (v1.0)",

      navbarMenu(
        "Genes",
        tabPanel(
          "Rankings",
          selectInput("select_gene_table_format", "Display:", 
                      choices = c('P-values', 'Ranks'), selected = 'Ranks'),
          withSpinner(dataTableOutput('mm25_genes'))
        ),
        tabPanel(
          "Visualize",
          fluidRow(
            column(
              width = 4,
              selectizeInput('select_gene', "Gene:", choices = NULL),
              helpText("Gene to visualize."),
              hr(),
              uiOutput('select_covariate'),
              helpText("Phenotype / covariate to compare gene expression or SNP counts against."),
              hr(),
              uiOutput('covariate_summary')
            ),
            column(
              width = 8,
              withSpinner(plotOutput('gene_plot', height = '800px'))
            )
          )
        )
      ),
      navbarMenu(
        "Pathways",
        tabPanel(
          "Ranking",
          selectInput("select_pathway_table_format", "Display:", 
                      choices = c('P-values', 'Ranks'), selected = 'Ranks'),
          withSpinner(dataTableOutput('mm25_pathways'))
        )
      ),
      navbarMenu(
        "Covariates",
        tabPanel(
          "P-value Distributions",
          selectInput('select_gene_pval_dist_type', "P-value type:",
                      choices = c("Adjusted (BH)", "Unadjusted"), selected = "Unadjusted"),
          withSpinner(plotOutput('gene_pval_dists', height = '4000px'))
        ),
        tabPanel(
          "Covariate Similarity",
          column(
            width = 3,
            selectInput('select_cov_similarity_feat_type', "Feature type:",
                        choices = c("Genes", "Pathways"), selected = "Genes"),
            selectInput('select_cov_similarity_cor_method', "Correlation type:",
                        choices = c("Pearson", "Spearman"), selected = "Pearson")
          ),
          column(
            width = 9,
            withSpinner(plotlyOutput('cov_similarity_heatmap', height = '800px'))
          )
        )
      ),
      # navbarMenu(
      #   "Functional Enrichment",
      #   tabPanel(
      #     "Phenotypes",
      #     #selectInput("select_fgsea_covariate", "Covariate:", choices =
               #     covariate_ids()),
      #     uiOutput("fgsea_select_covariate"),
      #     withSpinner(dataTableOutput("fgsea_results_indiv"))
      #   ),
      #   tabPanel(
      #     "MM25 Rankings",
      #     uiOutput("fgsea_select_ranking"),
      #     withSpinner(dataTableOutput("fgsea_results_combined"))
      #   ),
      #   tabPanel(
      #     "Summary",
      #     selectInput("select_fgsea_summary", "Display: ", choices = c('Covariates', 'Ranking Methods')),
      #     withSpinner(dataTableOutput("fgsea_summary_output"))
      #   )
      # )
      tabPanel(
        "Settings",
        selectInput("select_version", "Version:", choices=c('v1.0', 'v1.1'), selected = 'v1.1')
      )
    )
  )
}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
