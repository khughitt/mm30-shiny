#
# Feature weights shiny UI
# KH Jan 2020
#
library(DT)
library(feather)
library(ggdark)
library(ggforce)
library(gridExtra)
library(plotly)
library(shinythemes)
library(survival)
library(survminer)
library(tidyverse)
library(yaml)

options(stringsAsFactors = FALSE)
options(digits = 3)
set.seed(1)

# load MM25 combined scores
mm25_genes <- read_feather('/data/nih/fassoc/1.0/summary/gene_combined_scores.feather')
mm25_pathways <- read_feather('/data/nih/fassoc/1.0/summary/pathway_combined_scores.feather')                  

mm25_gene_pvals <- read_feather('/data/nih/fassoc/1.0/summary/gene_scores.feather')

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

# load individual feature / phenotype datasets (use non-redundant versions)
gene_infiles <- lapply(cfgs, function(x) { x$features$genes$rna })
gene_infiles <- sub('.feather', '_nr.feather', gene_infiles)

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

    dat <- tibble(
      phenotype = colnames(mm25_gene_pvals)[-1],
      pval = as.numeric(mm25_gene_pvals[mm25_gene_pvals$symbol == input$select_gene, -1])
    )

    dat <- dat %>%
      arrange(pval)
  })

  #
  # Form fields
  #
  output$select_covariate <- renderUI({
    req(input$select_gene)

    dat <- gene_pvals()

    # select choices (e.g. "MMRF IA14 overall survival (p = 0.003)")
    opts <- gene_pvals()$phenotype
    names(opts) <- sprintf("%s (p = %0.3f)",
                           gsub("_", " ", sub("_pval", "", dat$phenotype)),
                           dat$pval)

    selectInput('select_covariate', "Covariate:", choices = opts)
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


  #
  # Plots
  #
  output$gene_plot <- renderPlot({
    req(input$select_covariate)
    pheno <- input$select_covariate

    # work-around: detect MMRF prefix (includes underscore)
    # parts <- unlist(str_split(sub('_pval', '', pheno), '_'))
    pheno <- sub('_pval', '', pheno)

    if (startsWith(pheno, 'MMRF')) {
      dataset <- 'MMRF_IA14'
      covariate <- sub('MMRF_IA14_', '', pheno)

      # work-around 2...
      if (grepl('response', covariate)) {
        mdat <- 'treatment'
      } else if (grepl('survival', covariate)) {
        mdat <- 'survival' 
      } else {
        mdat <- 'clinical'
      }
    } else {
      dataset <- unlist(str_split(pheno, '_'))[1]
      covariate <- sub(paste0(dataset, '_'), '', pheno)

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

  output$gene_pval_dists <- renderPlot({
    set.seed(1)

    dat <- mm25_gene_pvals %>%
      sample_n(1000) %>%
      pivot_longer(-symbol, names_to = 'covariate', values_to ='pvalue')

    dat$dataset <- str_split(dat$covariate, '_', simplify = TRUE)[, 1]
    dat$dataset[dat$dataset == 'MMRF'] <- 'MMRF_IA14'
    dat$dataset <- factor(dat$dataset)

    dat$covariate <- factor(sub('_pval', '',  dat$covariate))

    dat$label <- trimws(gsub("_", " ", str_replace(dat$covariate, as.character(dat$dataset), '')))
    dat$label <- sprintf("%s (%s)", dat$label, dat$dataset)

    dat$method <- 'Logit'
    dat$method[grepl('survival', dat$covariate)] <- 'Survival'
    dat$method <- factor(dat$method)


    saveRDS(dat, "~/tmp.rds")

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
              plotOutput('gene_plot', height = '800px')
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
        "Datasets",
        tabPanel(
          "Gene P-values",
          plotOutput('gene_pval_dists', height = '4000px')
        )
      )
    )
  )
}




shinyApp(ui = ui, server = server, enableBookmarking = "url")
