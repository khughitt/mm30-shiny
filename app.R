#
# MM29 Shiny UI
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
library(survminer)
library(svglite)
library(tidyverse)
library(yaml)

source("R/plotting.R")

options(stringsAsFactors = FALSE)
options(spinner.color="#00bc8c")
options(scipen = 2, digits = 3)
set.seed(1)

# feature levels
feature_levels <- c(Genes = "genes", Pathways = "gene_sets")

# plot colors
#color_pal <- c("#36c764", "#c73699", "#3650c7", "#c7ac36", "#c7365f", "#bb36c7")
color_pal <- c("#00AFBB", "#E7B800", "#FC4E07","#BB3099","#EE0099","#0000AC")

# options for survival plot upper/lower expression cutoffs
surv_expr_cutoffs <- seq(5, 50, by = 5)
names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")

# load config
cfg <- read_yaml("config/config-v4.1.yml")

# data subset options ("all", "disease stage", etc.)
subset_opts <- names(cfg$mm25_scores$genes)
names(subset_opts) <- Hmisc::capitalize(gsub('_', ' ', subset_opts))

#############################
#
# Shiny Server
#
#############################
server <- function(input, output, session) {
  #
  # Reactives
  #

  # gene cytogenetic band mapping
  cyto_bands <- read_tsv(cfg$cytogenetic_bands, col_types = cols())

  # load dataset and covariate metadata
  phenotype_metadata <- reactive({
    message("[phenotype_metadata]")

    dataset_metadata <- read_tsv(cfg$geo_metadata, col_types = cols())

    infile <- cfg$phenotype_metadata
    message(sprintf("[phenotype_metadata] loading %s", infile))

    dat <- read_feather(infile)

    if ("feature_level" %in% colnames(dat)) {
      dat <- dat %>%
        filter(feature_level == "genes") %>%
        select(-feature_level)
    }

    # work-around: manually add MMRF to dataset metadata
    dataset_metadata <- rbind(dataset_metadata,
      c("MMRF", "MMRF CoMMpass Study IA18", NA, NA, NA, NA, NA, "GPL11154", NA,
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
    infile <- cfg$phenotype_clusters
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

    # add number of significant gene associations
    num_sig_assoc <- apply(mm25_feature_padjs(), 2, function (x) {
      sum(x < 0.01, na.rm = TRUE)
    })

    dat$num_sig_assoc <- num_sig_assoc[match(dat$covariate_id, names(num_sig_assoc))]

    dat
  })

  # load MM25 combined pvals
  mm25_gene_pvals_combined <- eventReactive(input$select_mm25_subset, {
    infile <- cfg$mm25_scores$genes[[input$select_mm25_subset]]

    message(sprintf("[mm25_gene_pvals_combined] loading %s", infile))

    dat <- read_feather(infile) %>%
      select(symbol, sumz_wt_pval, num_present, num_missing)

    dat$pval_adj <- p.adjust(dat$sumz_wt_pval, method='BH')

    dat %>%
      select(-sumz_wt_pval)
  })

  mm25_pathway_pvals_combined <- reactive({
    req(input$select_mm25_subset)

    infile <- cfg$mm25_scores$pathways[[input$select_mm25_subset]]

    message(sprintf("[mm25_pathway_pvals_combined] loading %s", infile))

    dat <- read_feather(infile)

    # add links to msigdb pathway info pages
    msigdb_ids <- sub("[^_]*_", "", dat$gene_set)

    msigdb_links <- sprintf("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/%s' target='_blank'>%s</a>",
                            msigdb_ids, dat$gene_set)
    dat <- dat %>%
      add_column(pathway = msigdb_links, .before = 1) %>%
      select(pathway, gene_set, sumz_wt_pval, num_present, num_missing)

    dat$pval_adj <- p.adjust(dat$sumz_wt_pval, method='BH')

    dat %>%
      select(-sumz_wt_pval)
  })

  # individual dataset p-values
  mm25_gene_pvals_indiv <- read_feather(cfg$gene_pvals_indiv)
  mm25_pathway_pvals_indiv <- read_feather(cfg$pathway_pvals_indiv)

  mm25_feature_pvals <- reactive({
    message("[mm25_feature_pvals]")

    if (rv$plot_feature_level == 'genes') {
      dat <- mm25_gene_pvals_indiv %>%
        rename(feature_id = symbol)
    } else {
      dat <- mm25_pathway_pvals_indiv %>%
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
    cfg_infiles <- Sys.glob(file.path(cfg$dataset_cfg_dir, "*.yml"))
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
    # mask <- startsWith(names(gene_infiles), "GSE")
    # gene_infiles[mask] <- sub(".feather", "_nr.feather", gene_infiles[mask])

    # work-around (may 31, 2021)
    # gene_infiles[mask] <- sub('/data/clean/geo/3.1', '/data/expr', gene_infiles[mask])
    # gene_infiles <- sub("/data/clean/mmrf/IA15/rnaseq/", "/data/expr/MMRF/", gene_infiles)

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

  # NOTE 2021-11-07:
  # disabling bulk load and switching to lazy-loading approach for now
  # to reduce memory requirements..
  # load individual dataset pathway-level expression data;
  # used to generate pathway-level feature vs. phenotype plots
  # pathway_data <- reactive({
  #   message("[pathway_data] init")
  #
  #   read_data <- function(x) {
  #     message(sprintf("[pathway_data] loading %s", x))
  #
  #     if (endsWith(x, 'parquet')) {
  #       read_parquet(x)
  #     } else {
  #       read_feather(x)
  #     }
  #   }
  #
  #   infiles <- lapply(dataset_cfgs(), function(x) {
  #     x$features$gene_sets$rna
  #   })
  #
  #   dat <- lapply(infiles, read_data)
  #   names(dat) <- names(dataset_cfgs())
  #
  #   dat
  # })

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


  feature_plot_dat <- reactive({
    req(rv$plot_covariate)
    req(input$select_survival_expr_cutoffs)

    message("[feature_plot_dat]")

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
      # load pathway-level expression data
      infile <- dataset_cfgs()[[dataset_id]]$features$gene_sets$rna

      if (endsWith(infile, 'parquet')) {
        pathway_expr <- read_parquet(infile)
      } else {
        pathway_expr <- read_feather(infile)
      }

      #feat_dat <- pathway_data()[[dataset_id]] %>%
      feat_dat <- pathway_expr %>%
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

      data.frame(feature = feat_dat, time = time_dat, event = event_dat)
    } else if (assoc_method == "logit") {
      # violin plot
      response <- factor(pull(pheno_data()[[dataset_id]], assoc_cfg$params$field))
      data.frame(feature = feat_dat, response)
    } else if (assoc_method == "deseq") {
      # TODO/NEXT STEPS...
      # (DESeq plot?..)
      data.frame()
    }
  })

  #
  # html output
  #
  output$covariate_summary <- renderUI({
    pheno_mdata <- selected_phenotype_metadata()

    tag_list <- tagList(
      tags$b(pheno_mdata$title),
      br(),
      tags$b(tags$a(href=pheno_mdata$urls, target='_blank', pheno_mdata$dataset)),
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

    names(opts) <- sprintf("%s (padj = %0.3f)",
                           gsub("_", " ", sub("_pval", "", dat$phenotype)),
                           dat$padj)

    selectInput("select_plot_covariate", "Covariate:", choices = opts)
  })


  rv <- reactiveValues(
    # version_subset = NULL,
    plot_feature_level = NULL,
    plot_feature = NULL,
    plot_covariate = NULL
  )

  # manually trigger subset selection
  # updateSelectInput(session, "select_mm25_subset", "Data Subset:",
  #                   choices = subset_opts, selected = "all")

  # Using static set of GRCh38 genes from annotables; the datasets in MM25 also
  # include some other gene symbols, but there are generally quite rare and not
  # likely to be of interest
  updateSelectizeInput(session, "select_coex_gene1", choices = grch38$symbol,
                       selected = "MCL1", server = TRUE)

  updateSelectizeInput(session, "select_coex_gene2", choices = grch38$symbol,
                       selected = "PBXIP1", server = TRUE)

  observe({
    updateSelectInput(session, "select_coex_experiment", choices=names(dataset_cfgs()))
  })

  #
  # Event Hanlders
  #

  # version subset
  observeEvent(input$select_mm25_subset, {
    if (input$select_mm25_subset != 'loading...') {
      rv$plot_feature_level <- "genes"

      message("observeEvent(input$select_mm25_subset)")

      updateSelectInput(session, "select_plot_feature_level", "Feature Level:",
                        choices = feature_levels, selected = "genes")
    }
  })

  observeEvent(input$select_plot_feature_level, {
    req(input$select_mm25_subset)

    # save selected dataset/covariate
    prev_covariate = rv$plot_covariate

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

    # set initial covariate & feature to trigger plot update (still not behaving as
    # expected.. have to change covariate after switching to pathways for plot to be
    # rendered..)
    rv$plot_covariate <- prev_covariate
    rv$plot_feature = select_choices[1]
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
  })

  #
  # Tables
  #
  output$mm25_gene_pvals_combined_table <- renderDataTable({
    req(input$select_mm25_subset)

    # float_cols <- paste0(c("min", "sumlog", "sumz", "sumz_wt"), "_pval")
    float_cols <- c("pval_adj")

    dat <- mm25_gene_pvals_combined()

    # add biotype and description
    gene_annot <- grch38 %>%
      select(symbol, description, biotype) %>%
      group_by(symbol) %>%
      slice(1) %>%
      ungroup()

    # remove the source info from gene description
    gene_annot$description <- str_split(gene_annot$description, " \\[", simplify = TRUE)[, 1]

    dat <- dat %>%
      left_join(gene_annot, by = 'symbol') %>%
      select(symbol, biotype, everything())

    # add gene cytogenetic bands
    dat <- dat %>%
      inner_join(cyto_bands, by = 'symbol')

    # construct data table
    out <- DT::datatable(dat, style = "bootstrap", options = list(pageLength = 15))

    out <- out %>%
        formatSignif(columns = float_cols, digits = 3)

  })

  output$mm25_pathway_pvals_combined_table <- renderDataTable({
    req(input$select_mm25_subset)
    # req(input$select_pathway_table_format)

    message("mm25_pathway_pvals_combined_table")

    # float_cols <- paste0(c("min", "sumlog", "sumz", "sumz_wt"), "_pval")
    float_cols <- c("pval_adj")

    dat <- mm25_pathway_pvals_combined()

    # construct data table
    DT::datatable(dat %>% select(-gene_set),
                         style = "bootstrap", escape = FALSE,
                         options = list(pageLength = 15)) %>%
      formatSignif(columns = float_cols, digits = 3)
  })

  output$download_plot_data <- downloadHandler(
    filename = function() {
      # filename parts
      pheno <- rv$plot_covariate
      pheno <- sub("_pval", "", pheno)

      sprintf("%s.tsv", pheno)
    },
    content = function(file) {
      write_tsv(feature_plot_dat(), file)
    }
  )

  output$download_plot_svg <- downloadHandler(
      filename = function() {
        # filename parts
        pheno <- rv$plot_covariate
        pheno <- sub("_pval", "", pheno)

        sprintf("%s.svg", pheno)
      },

      content = function(file) {
        svglite(file)
        print(getFeatPlot())
        dev.off()
      }
  )

  output$download_plot_tiff <- downloadHandler(
      filename = function() {
        # filename parts
        pheno <- rv$plot_covariate
        pheno <- sub("_pval", "", pheno)

        sprintf("%s.tiff", pheno)
      },

      content = function(file) {
        tiff(file, width=6, height=5, res=300, units="in")
        print(getFeatPlot())
        dev.off()
      }
  )

  #
  # Plots
  #
  featurePlot <- reactive({
    req(feature_plot_dat)
    req(rv$plot_covariate)

    message("[feature_plot]")

    pheno <- rv$plot_covariate
    pheno <- sub("_pval", "", pheno)

    # get dataset and covariate names
    dataset_id <- unlist(str_split(pheno, "_"))[1]
    covariate <- sub(paste0(dataset_id, "_"), "", pheno)

    # get config for selected feature-phenotype association
    assoc_cfg <- dataset_cfgs()[[dataset_id]]$phenotypes$associations[[covariate]]
    assoc_method <- assoc_cfg$method

    dat <- feature_plot_dat()

    if (assoc_method == "survival") {
      # survival plot
      plot_survival(dat, dataset_id, covariate, rv$plot_feature,
                    input$select_survival_expr_cutoffs, color_pal)
    } else if (assoc_method == "logit") {
      # violin plot
      plot_categorical(dat, dataset_id, covariate, rv$plot_feature, color_pal)
    } else if (assoc_method == "deseq") {
      # TODO...
      #plot_deseq(dat, dataset_id, covariate, rv$plot_feature, color_pal)
    }
  })
  output$feature_plot <- renderPlot(featurePlot())

  getFeatPlot <- function() {
    featurePlot()
  }

  output$gene_coex_plot <- renderPlot({
    req(gene_data())
    req(input$select_coex_gene1)
    req(input$select_coex_gene2)

    gene1 <- input$select_coex_gene1
    gene2 <- input$select_coex_gene2

    dat_name <- input$select_coex_experiment

    dat <- gene_data()[[dat_name]] %>%
      filter(symbol %in% c(gene1, gene2))

    if (nrow(dat) == 2) {
      dat_long <- dat %>%
        column_to_rownames('symbol') %>%
        t() %>%
        as.data.frame()

      # compute pearson correlation
      coex_cor <- cor(dat_long[, 1], dat_long[, 2], use = 'pairwise.complete')

      ggplot(dat_long, aes_string(gene1, gene2)) +
          geom_point() +
          geom_smooth(method = lm) +
          ggtitle(sprintf("%s (cor = %0.2f)", dat_name, coex_cor))
    } else {
      # if both genes cannot be found in the dataset, for now, don't do anything..
      NULL
    }
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
      theme = shinytheme("flatly"),
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
            downloadButton("download_plot_svg", "Download (.svg)"),
            downloadButton("download_plot_tif", "Download (.tiff)"),
            downloadButton("download_plot_data", "Download (.tsv)"),
            withSpinner(plotOutput("feature_plot", height = "740px"))
          )
        )
      ),
      tabPanel(
        "Co-expression",
        column(
          width = 3,
          selectizeInput("select_coex_gene1", "Gene 1", choices = NULL),
          selectizeInput("select_coex_gene2", "Gene 2", choices = NULL),
          selectizeInput("select_coex_experiment", "Experiment", choices = NULL)

        ),
        column(
          width = 9,
          withSpinner(plotOutput("gene_coex_plot", height = "740px"))
        )
      ),
      tabPanel(
        "Settings",
        selectInput("select_mm25_subset", "Subset:", 
                    choices = subset_opts, selected = "all")
      )
    )
  )
}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
