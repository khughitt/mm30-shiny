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
library(survival)
library(survminer)
library(tidyverse)
library(yaml)

options(stringsAsFactors = FALSE)
options(digits = 3)
set.seed(1)

# load MM25 combined scores / pvals

# genes
mm25_gene_scores <- read_feather('/data/nih/fassoc/1.0/summary/gene_scores.feather')

mm25_genes <- read_feather('/data/nih/fassoc/1.0/summary/gene_pvals.feather') %>%
  inner_join(mm25_gene_scores, by = 'symbol')

# pathways
mm25_pathway_scores <- read_feather('/data/nih/fassoc/1.0/summary/pathway_scores.feather')

mm25_pathways <- read_feather('/data/nih/fassoc/1.0/summary/pathway_pvals.feather') %>%
  inner_join(mm25_pathway_scores, by = 'gene_set')

# clean-up
rm(mm25_gene_scores, mm25_pathway_scores)

# individual dataset p-values
mm25_gene_pvals <- read_feather('/data/nih/fassoc/1.0/summary/gene_pvals_indiv.feather')

mm25_gene_padjs <- mm25_gene_pvals %>%
    mutate_at(vars(-symbol), p.adjust, method = 'BH')

# list of covariates
covariate_ids <- colnames(mm25_gene_pvals)[-1]  

# add links to msigdb pathway info pages
msigdb_ids <- sub("[^_]*_", "", mm25_pathways$gene_set) 

msigdb_links <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>", 
                        msigdb_ids, mm25_pathways$gene_set)
mm25_pathways <- mm25_pathways %>%
  add_column(pathway = msigdb_links, .before = 1)

# load dataset configs
cfg_infiles <- Sys.glob('~/d/r/nih/fassoc/datasets/*.yml')
cfgs <- lapply(cfg_infiles, read_yaml)
names(cfgs) <- toupper(sub('.yml', '', basename(cfg_infiles)))

# load individual feature / phenotype datasets
gene_infiles <- lapply(cfgs, function(x) { x$features$genes$rna })

# for GEO datasets, make sure to use non-redundant (nr) versions of expression data
mask <- startsWith(names(gene_infiles), "GSE")
gene_infiles[mask] <- sub('.feather', '_nr.feather', gene_infiles[mask])

gene_data <- lapply(gene_infiles, read_feather)
names(gene_data) <- names(cfgs)

pheno_data <- list()

for (dataset in names(cfgs)) {
  cfg <- cfgs[[dataset]]

  pheno_data[[dataset]] <- list()

  for (phenotype in names(cfg$phenotypes)) {
    infile <- cfg$phenotypes[[phenotype]]$path

    if (endsWith(infile, 'tsv') || endsWith(infile, 'tsv.gz')) {
      pheno_data[[dataset]][[phenotype]] <- read_tsv(infile, col_types = cols())
    } else if (endsWith(infile, 'feather')) {
      pheno_data[[dataset]][[phenotype]] <- read_feather(infile)
    }
  }
}

# columns to round
float_cols <- colnames(mm25_genes)[!colnames(mm25_genes) %in%
                                   c('symbol', 'num_present', 'num_missing')]

# ggplot theme
theme_dark <- dark_theme_gray() +
  theme(axis.text.x = element_text(angle = 90),
        plot.background = element_rect(fill = "#222222"),
        panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
        panel.grid.major = element_line(color = "#555555", size = 0.2),
        panel.grid.minor = element_line(color = "#555555", size = 0.2))


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
  gene_pvals <- reactive({
    req(input$select_gene)

    mask <- mm25_gene_pvals$symbol == input$select_gene

    dat <- tibble(
      phenotype = covariate_ids,
      pval = as.numeric(mm25_gene_pvals[mask, -1]),
      padj = as.numeric(mm25_gene_padjs[mask, -1])
    )

    dat <- dat %>%
      arrange(pval)
  })

  fgsea_summary_indiv <- reactive({
    res <- read_tsv("/data/nih/fgsea/1.0-orig/results/summary/mm25_pvals_indiv_summary.tsv", 
                    col_types = cols())
    res %>%
      arrange(desc(num_sig))
  })

  fgsea_summary_combined <- reactive({
    pvals <- read_tsv("/data/nih/fgsea/1.0-orig/results/summary/mm25_pvals_combined_summary.tsv", 
                      col_types = cols())

    scores <- read_tsv("/data/nih/fgsea/1.0-orig/results/summary/mm25_combined_summary.tsv",
                       col_types = cols())

    rbind(pvals, scores) %>%
      arrange(desc(num_sig))
  })

  fgsea_results_indiv <- reactive({
    read_parquet("/data/nih/fgsea/1.0-orig/results/merged/mm25_pvals_indiv.parquet") %>%
        select(-dataset) %>%
        filter(padj < 0.05) %>%
        arrange(padj)
  })

  fgsea_results_indiv_filtered <- reactive({
    req(input$select_fgsea_covariate)

    dat <- fgsea_results_indiv() %>%
        filter(field == input$select_fgsea_covariate) %>%
        select(-field)

    dat$pathway <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>", 
                            dat$pathway, dat$pathway)

    dat
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

  output$fgsea_select_covariate <- renderUI({
    fgsea_summary <- fgsea_summary_indiv()
    
    opts <- covariate_ids

    num_sig <- fgsea_summary$num_sig[match(opts, fgsea_summary$field)]

    # include number of significant terms in labels
    names(opts) <- sprintf("%s (#sig: %d)", opts, num_sig)

    # order labels in decreasing order of functional enrichment
    opts <- opts[match(fgsea_summary$field, opts)]

    selectInput("select_fgsea_covariate", "Covariate:", choices = opts)
  })

  updateSelectizeInput(session, 'select_gene', 
                       choices = as.character(mm25_genes$symbol), server = TRUE)

  # 
  # Tables
  #
  output$mm25_genes <- renderDataTable({
    DT::datatable(mm25_genes, style = 'bootstrap', options = list(pageLength = 15)) %>%
      formatRound(columns = float_cols, digits = 5)
  })
  output$mm25_pathways <- renderDataTable({
    DT::datatable(mm25_pathways %>% select(-gene_set), 
                  style = 'bootstrap', escape = FALSE, options = list(pageLength = 15)) %>%
      formatRound(columns = float_cols, digits = 5)
  })

  output$fgsea_results_indiv <- renderDataTable({
    DT::datatable(fgsea_results_indiv_filtered(), 
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

    if (startsWith(pheno, 'MMRF')) {
      # work-around...
      if (grepl('response', covariate)) {
        mdat <- 'treatment'
      } else if (grepl('survival', covariate)) {
        mdat <- 'survival' 
      } else {
        mdat <- 'clinical'
      }
    } else {
      mdat <- names(cfgs[[dataset]]$phenotypes)[1]
    }

    # config
    cfg <- cfgs[[dataset]]$phenotypes[[mdat]]$associations[[covariate]]

    if (cfg$method == 'survival') {
      # survival data
      feature <- gene_data[[dataset]] %>%
        filter(symbol == input$select_gene) %>%
        select(-symbol) %>%
        as.numeric()

      time_dat <- pull(pheno_data[[dataset]][[mdat]], cfg$params$time)
      event_dat <- pull(pheno_data[[dataset]][[mdat]], cfg$params$event)

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

      # perform fit on binarized data
      fit <- survfit(Surv(time, event) ~ Expression, data = dat)

      # display a kaplan meier plot for result
      ggsurvplot(fit, data = dat, ggtheme = theme_dark, palette = color_pal)
    } else {
      # logistic regression fit plot
      feature <- gene_data[[dataset]] %>%
        filter(symbol == input$select_gene) %>%
        select(-symbol) %>%
        as.numeric()

      # pheno data
      response <- pull(pheno_data[[dataset]][[mdat]], cfg$params$field)

      dat <- data.frame(feature, response)

      dat$response <- factor(dat$response)

      # drop any entries with missing values
      dat <- dat[!is.na(dat$response), ]

      # draw violin + jitter plot 
      set.seed(1)

      ggplot(dat, aes(x = response, y = feature)) +
        # geom_boxplot(aes(fill = response, color = response), alpha = 0.7, outlier.shape = NA) +
        geom_violin(aes(fill = response, color = response), alpha = 0.8, draw_quantiles = c(0.5)) +
        geom_jitter(aes(color = response), alpha = 0.9) +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        ggtitle(sprintf("%s: %s vs. %s", dataset, input$select_gene, covariate)) +
        xlab(covariate) +
        ylab(sprintf("%s expression", input$select_gene)) +
        dark_theme_gray() +
        theme(axis.text.x = element_text(angle = 90),
              plot.background = element_rect(fill = "#222222"),
              panel.border = element_rect(colour = "#333333", fill = NA, size = 1),
              panel.grid.major = element_line(color = "#555555", size = 0.2),
              panel.grid.minor = element_line(color = "#555555", size = 0.2))
    }
  })

  output$covariate_similarity <- renderPlotly({
    # determine covariate functional groups from labels

    cnames <- colnames(mm25_gene_pvals)[-1]

    cov_labels <- rep("", ncol(mm25_gene_pvals) - 1)

    cov_labels[grepl('survival|died', cnames)] <- 'survival'
    cov_labels[grepl('treatment|response', cnames)] <- 'treatment'
    cov_labels[grepl('status|stage|ecog|pfs_event|relapsed', cnames)] <- 'stage'

    cov_labels <- factor(cov_labels)

    annots <- data.frame(type = cov_labels)
    ann_colors <- list(type = c('#799d15', '#15799e', '#9e1579'))

    # generate covariate correlation matrix
    cor_mat <- cor(mm25_gene_pvals[, -1], method = 'pearson', 
                   use = 'pairwise.complete.obs')

    pal <- color_pal[1:3]
    names(pal) <- c('survival', 'treatment', 'stage')

    heatmaply(cor_mat, row_side_colors = annots, row_side_palette = pal)
  })

  output$gene_pval_dists <- renderPlot({
    req(input$select_gene_pval_dist_type)

    set.seed(1)

    # adjusted / unadjusted p-values
    if (input$select_gene_pval_dist_type == "Adjusted (BH)") {
        dat <- mm25_gene_padjs  %>%
        sample_n(1000) %>%
        pivot_longer(-symbol, names_to = 'covariate', values_to ='pvalue')
    } else {
        dat <- mm25_gene_pvals %>%
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


    npages <- ceiling((ncol(mm25_gene_pvals) - 1) / 9)

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
        HTML("#title { font-family: 'Roboto Mono', monospace; }")
      )
    ),

    navbarPage(
      id = "main",
      theme = shinytheme("darkly"),
      title = "MM25",
      windowTitle = "MM25",

      navbarMenu(
        "Genes",
        tabPanel(
          "Rankings",
          dataTableOutput('mm25_genes')
        ),
        tabPanel(
          "Visualize",
          fluidRow(
            column(
              width = 3,
              selectizeInput('select_gene', "Gene:", choices = NULL),
              uiOutput('select_covariate')
            ),
            column(
              width = 9,
              plotOutput('gene_plot', height = '960px')
            )
          )
        )
      ),
      navbarMenu(
        "Pathways",
        tabPanel(
          "Ranking",
          dataTableOutput('mm25_pathways')
        )
      ),
      navbarMenu(
        "Covariates",
        tabPanel(
          "P-value Distributions",
          selectInput('select_gene_pval_dist_type', "P-value type:", 
                      choices = c("Adjusted (BH)", "Unadjusted"), selected = "Unadjusted"),
          plotOutput('gene_pval_dists', height = '4000px')
        ),
        tabPanel(
          "Covariate Similarity",
          plotlyOutput('covariate_similarity', height = '800px')
        )
      ),
      navbarMenu(
        "Functional Enrichment",
        tabPanel(
          "Covariates",
          #selectInput("select_fgsea_covariate", "Covariate:", choices = covariate_ids),
          uiOutput("fgsea_select_covariate"),
          dataTableOutput("fgsea_results_indiv")
        ),
        tabPanel(
          "Summary",
          selectInput("select_fgsea_summary", "Display: ", choices = c('Covariates', 'Ranking Methods')),
          dataTableOutput("fgsea_summary_output")
        )
      )
    )
  )
}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
