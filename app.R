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
library(survminer)
library(svglite)
library(tidyverse)
library(yaml)

source("R/plotting.R")

options(stringsAsFactors = FALSE)
options(spinner.color = "#00bc8c")
options(scipen = 2, digits = 3)
set.seed(1)

# TESTING (Sep23)
data_dir <- "../data"

# feature levels
feature_types <- c(Genes = "genes", Pathways = "pathways")

# load config
cfg <- read_yaml("config/config-v7.0.yml")

# table options
tableOpts <- list(pageLength = 15)

# human-readable column names (includes placeholders for hidden names);
# hidden fields are used to store unformatted values
geneColNames <- c("Symbol", "Biotype", "Chr", "Cell Cycle Phase",
                  "HIDDEN1", "HIDDEN2", "HIDDEN3", "HIDDEN4", "HIDDEN5",
                  "Total", "Disease Stage", "Survival (OS)",
                  "Survival (PFS)", "Treatment Response")

covariateOpts <- c("All" = "all", "Disease Stage" = "disease_stage", "Survival (OS)" = "survival_os",
                   "Survival (PFS)" = "survival_pfs",
                   "Treatment Response" = "treatment_response")

# options for survival plot upper/lower expression cutoffs
surv_expr_cutoffs <- cfg$surv_cutoff_opts
names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")

# load mm30 gene/pathways scores
mm30_genes <- read_feather(file.path(data_dir, "mm30", "mm30_gene_scores.feather"))

mm30_pathways <- read_feather(file.path(data_dir, "mm30", "mm30_pathway_scores.feather"))

mm30_pathways$pathway_html <- sprintf("<a href='%s?geneSetName=%s' target='_blank'>%s</a>",
                                      "http://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp",
                                      mm30_pathways$pathway, mm30_pathways$pathway)

# load metadata
gene_mdata <- read_tsv(file.path(data_dir, "metadata/genes.tsv"), col_types = cols()) %>%
  select(-chr_subband)

covariates <- read_yaml(file.path(data_dir, "metadata/covariates.yml"))

# load mm30 dataset column metadata
col_mdata <- list()

for (acc in names(covariates)) {
  if (acc == "MMRF") {
    infile <- file.path(data_dir, "mmrf/column-metadata.feather")
  } else {
    infile <- file.path(data_dir, sprintf("geo/%s/column-metadata.feather", acc))
  }
  col_mdata[[acc]] <- read_feather(infile)
}

# load mm30 individual GEO/MMRF gene expression datasets
geo_dir <- file.path(data_dir, "geo")

mm30_expr <- list()

for (acc in list.files(geo_dir)) {
  mm30_expr[[acc]] <- read_feather(file.path(geo_dir, acc, "data.feather"))
}
mm30_expr[["MMRF"]] <- read_feather(file.path(data_dir, "mmrf", "data.feather"))

# load dataset + covariate metadata
covariate_mdata <- read_feather(file.path(data_dir, "metadata/covariates.feather"))

# individual dataset p-values
mm30_gene_pvals_indiv <- read_feather(file.path(data_dir,
                                                "fassoc/gene_association_pvals.feather"))
mm30_pathway_pvals_indiv <- read_feather(file.path(data_dir,
                                                "fassoc/pathway_association_pvals.feather"))

# table headers (colspan include 5 hidden columns)
geneTableHeader <- htmltools::withTags(table(
  class = "display",
  thead(
    tr(
      th(rowspan = 2, "Symbol"),
      th(rowspan = 2, "Biotype"),
      th(rowspan = 2, "Chr"),
      th(rowspan = 2, "Cell Cycle Phase"),
      th(colspan = 8, "Gene Associations (Rank / Padj)"),
    ),
    tr(
       th("All Covariates"),
       th("Disease Stage"),
       th("Survival (OS)"),
       th("Survival (PFS)"),
       th("Treatment Response")
    )
  )
))

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
  # Reactives
  #
  gene_df <- reactive({
    # add biotype and description
    gene_annot <- grch38 %>%
      select(symbol, description, biotype) %>%
      group_by(symbol) %>%
      slice(1) %>%
      ungroup()

    # remove the source info from gene description
    gene_annot$description <- str_split(gene_annot$description, " \\[", simplify = TRUE)[, 1]

    dat <- mm30_genes %>%
      select(-num_datasets) %>%
      left_join(gene_annot, by = "symbol")

    # add gene annotations
    dat <- dat %>%
      left_join(gene_mdata, by = "symbol")

    # add covariate ranks
    dat$all_rank <- rank(dat$all, ties.method = "first")
    dat$disease_stage_rank <- rank(dat$disease_stage, ties.method = "first")
    dat$survival_os_rank <- rank(dat$survival_os, ties.method = "first")
    dat$survival_pfs_rank <- rank(dat$survival_pfs, ties.method = "first")
    dat$treatment_rank <- rank(dat$treatment_response, ties.method = "first")

    dat <- dat %>%
      select(-description) %>%
      select(symbol, biotype, chr_region, cell_cycle_phase, everything())

    dat$all_formatted <- sprintf("%-6d (%s)", dat$all_rank, format(dat$all))

    dat$disease_stage_formatted <- sprintf("%-6d (%s)",
                                           dat$disease_stage_rank,
                                           format(dat$disease_stage))

    dat$survival_os_formatted <- sprintf("%-6d (%s)",
                                      dat$survival_os_rank,
                                      format(dat$survival_os))

    dat$survival_pfs_formatted <- sprintf("%-6d (%s)",
                                      dat$survival_pfs_rank,
                                      format(dat$survival_pfs))


    dat$treatment_formatted <- sprintf("%-6d (%s)", dat$treatment_rank, format(dat$treatment_response))

    dat %>% select(-disease_stage_coef, -survival_coef, -all_rank, -survival_os_rank,
                   -survival_pfs_rank, -disease_stage_rank, -treatment_rank)
  })

  mm30_plot_pvals <- reactive({
    message("[mm30_plot_pvals]")

    if (rv$plot_feature_type == "genes") {
      dat <- mm30_gene_pvals_indiv %>%
        rename(feature_id = symbol)
    } else {
      dat <- mm30_pathway_pvals_indiv %>%
        rename(feature_id = pathway)
    }
    dat
  })

  mm30_plot_padjs <- reactive({
    message("[mm30_plot_padjs]")

    mm30_plot_pvals() %>%
      mutate_at(vars(-feature_id), p.adjust, method = "BH")
  })

  # load gene-specific co-ex data
  gene_coex <- reactive({
    req(input$coex_gene1)
    gene <- input$coex_gene1
    expt_id <- input$coex_experiment

    message("[gene_coex]")

    infile <- file.path(data_dir, sprintf("coex/%s.feather", gene))

    if (!file.exists(infile)) {
      return(data.frame())
    }

    read_feather(infile) %>%
      filter(dataset == expt_id) %>%
      select(-dataset)
  })

  selected_feature_pvals <- reactive({
    req(rv$plot_feature)

    message("[selected_feature_pvals]")

    mask <- mm30_plot_pvals()$feature_id == rv$plot_feature

    covariate_ids <- colnames(mm30_plot_pvals())[-1]

    dat <- tibble(
      phenotype = covariate_ids,
      pval = as.numeric(mm30_plot_pvals()[mask, -1]),
      padj = as.numeric(mm30_plot_padjs()[mask, -1])
    )

    dat %>%
      arrange(pval)
  })

  selected_covariate_mdata <- reactive({
    req(rv$plot_covariate)

    message("[selected_covariate_mdata]")

    # determine selected dataset / covariate
    dataset_id <- unlist(str_split(rv$plot_covariate, "_"))[1]
    pheno <- sub("_pval", "", sub(paste0(dataset_id, "_"), "", rv$plot_covariate))

    # retrieve relevant metadata
    covariate_mdata %>%
      filter(dataset ==  dataset_id & phenotype == pheno)
  })

  feature_plot_dat <- reactive({
    req(rv$plot_covariate)
    req(input$survival_cutoff_percentile)

    message("[feature_plot_dat]")

    pheno <- rv$plot_covariate

    # get dataset and covariate names
    dataset_id <- unlist(str_split(pheno, "_"))[1]
    covariate <- sub(paste0(dataset_id, "_"), "", pheno)

    # feature data
    if (rv$plot_feature_type == "genes") {
      feat_dat <- mm30_expr[[dataset_id]] %>%
        filter(symbol == rv$plot_feature) %>%
        select(-symbol) %>%
        as.numeric()
    } else {
      # load pathway-level expression data
      message("......")
      message(dataset_id)
      message("......")

      data_subdir <- ifelse(dataset_id == "MMRF", "mmrf", file.path("geo", dataset_id))

      infile <- file.path(data_dir, data_subdir, "data_pathways.feather")
      pathway_expr <- read_feather(infile)

      #feat_dat <- pathway_data()[[dataset_id]] %>%
      feat_dat <- pathway_expr %>%
        filter(pathway == rv$plot_feature) %>%
        select(-pathway) %>%
        as.numeric()
    }

    # get config for selected feature-phenotype association
    covariate_info <- covariates[[dataset_id]][[covariate]]

    if (is.null(covariate_info)) {
      stop("Error retrieving covariate info!")
    }

    assoc_method <- covariate_info$method

    if (assoc_method == "survival") {
      # survival plot
      pheno_dat <- col_mdata[[dataset_id]]

      time_dat <- pull(pheno_dat, covariate_info$params$time)
      event_dat <- pull(pheno_dat, covariate_info$params$event)

      data.frame(feature = feat_dat, time = time_dat, event = event_dat)
    } else if (assoc_method == "logit") {
      # violin plot
      response <- factor(pull(col_mdata[[dataset_id]], covariate_info$params$field))
      data.frame(feature = feat_dat, response)
    } else if (assoc_method == "deseq") {
      # violin plot?
      model_parts <- unlist(strsplit(covariate_info$design$full, " "))

      var1 <- substring(model_parts[[1]], 2)
      feat1 <- factor(pull(col_mdata[[dataset_id]], var1))

      # ex. ~treatment + replicate
      if (length(model_parts) == 3) {
        dat <- data.frame(feature = feat_dat, feat1)
        colnames(dat)[2] <- var1
      } else if (length(model_parts) == 5) {
        # ex. ~time_hours + treatment + replicate
        var2 <- model_parts[[3]]
        feat2 <- factor(pull(col_mdata[[dataset_id]], var2))
        dat <- data.frame(feature = feat_dat, feat1, feat2)
        colnames(dat)[2:3] <- c(var1, var2)
      }
      dat
    }
  })

  #
  # html output
  #
  output$gene_info <- renderUI({
    message("output$gene_info")

    tag_list <- tagList(
      tags$b(rv$selected_gene),
      br()
    )

    tags$div(tag_list)
  })

  output$covariate_summary <- renderUI({
    message("output$covariate_summary")

    mdata <- selected_covariate_mdata()

    tag_list <- tagList(
      tags$b(mdata$title),
      br(),
      tags$b(tags$a(href = mdata$urls, target = "_blank", mdata$dataset)),
      br()
    )

    if (!is.na(mdata$name)) {
      authors <- str_to_title(str_match(mdata$name, "[[:alpha:]]+"))
      year <- str_split(mdata$submission_date, " ", simplify = TRUE)[, 3]

      attribution <- sprintf("%s (%s)", authors, year)
      tag_list <- tagAppendChildren(tag_list, attribution, br())
    }

    tag_list <- tagAppendChildren(
      tag_list,
      "# Samples:",
      tags$b(mdata$num_samples),
      br(),
      "# Genes:",
      tags$b(mdata$num_genes),
      br()
    )

    if (!is.na(mdata$overall_design)) {
      tag_list <- tagAppendChildren(tag_list, hr(),
                                    tags$b("Overall Design:"), br(),
                                    mdata$abstract, br())
    }

    tags$div(tag_list)
  })

  #
  # Form fields
  #
  output$plot_covariate <- renderUI({
    req(rv$plot_feature)

    message("output$plot_covariate")

    dat <- selected_feature_pvals()

    # select choices (e.g. "MMRF overall survival (p = 0.003)")
    opts <- selected_feature_pvals()$phenotype

    names(opts) <- sprintf("%s (padj = %0.3f)",
                           gsub("_", " ", sub("_pval", "", dat$phenotype)),
                           dat$padj)
    # TODO (Oct 9, 2023): Add dropdown menu to choose subset for plots...

    selectInput("plot_covariate", "Covariate:", choices = opts)
  })


  rv <- reactiveValues(
    plot_feature_type = NULL,
    plot_feature = NULL,
    plot_covariate = NULL,
    selected_gene = "",
  )

  # manually trigger subset selection
  # updateSelectInput(session, "select_mm30_subset", "Data Subset:",
  #                   choices = subset_opts, selected = "all")

  # Using static set of GRCh38 genes from annotables; the datasets in mm30 also
  # include some other gene symbols, but there are generally quite rare and not
  # likely to be of interest
  updateSelectizeInput(session, "coex_gene1", choices = grch38$symbol,
                       selected = "MCL1", server = TRUE)

  updateSelectizeInput(session, "coex_gene2", choices = grch38$symbol,
                       selected = "PBXIP1", server = TRUE)

  observe({
    updateSelectInput(session, "coex_experiment",  choices = names(covariates))
  })

  #
  # Event Hanlders
  #
  observeEvent(input$mm30_gene_pvals_combined_table_row_last_clicked, {
    ind <- input$mm30_gene_pvals_combined_table_row_last_clicked
    target_gene <- gene_df()$symbol[ind]

    rv$selected_gene <- ifelse(rv$selected_gene == target_gene, "", target_gene)
  })


  observeEvent(input$feature_type, {
    message("observeEvent(input$feature_type)")

    # save selected dataset/covariate
    prev_covariate <- rv$plot_covariate

    rv$plot_feature <- NULL
    rv$plot_covariate <- NULL

    message(input$feature_type)

    if (input$feature_type != "") {
      rv$plot_feature_type <- input$feature_type
    }

    if (rv$plot_feature_type == "genes") {
      select_choices <- as.character(mm30_genes$symbol)
    } else {
      select_choices <- as.character(mm30_pathways$pathway)
    }

    updateSelectizeInput(session, "feature", choices = select_choices, server = TRUE)

    # set initial covariate & feature to trigger plot update (still not behaving as
    # expected.. have to change covariate after switching to pathways for plot to be
    # rendered..)
    rv$plot_covariate <- prev_covariate
    rv$plot_feature <- select_choices[1]
  })

  observeEvent(input$plot_covariate, {
    message("observeEvent(input$plot_covariate)")
    rv$plot_covariate <- input$plot_covariate
  })

  observeEvent(input$feature, {
    req(input$feature_type)

    message("observeEvent(input$feature)")

    rv$plot_covariate <- NULL
    rv$plot_feature <- input$feature
  })

  #
  # Tables
  #
  output$mm30_gene_pvals_combined_table <- renderDataTable({
    req(gene_df())

    message("output$mm30_gene_pvals_combined_table")

    opts <- list(
      pageLength = 15,
      columnDefs = list(
        list(orderData = 4, targets = 9),
        list(orderData = 5, targets = 10),
        list(orderData = 6, targets = 11),
        list(orderData = 7, targets = 12),
        list(visible = FALSE, targets = 4:8)
      )
    )

    DT::datatable(gene_df(),
                  style = "bootstrap",
                  escape = FALSE,
                  colnames = geneColNames,
                  rownames = FALSE,
                  selection = "single",
                  options = opts)

    # disabling for now; conflicts with columnDefs used for sorting formatting fields
    #container = geneTableHeader,
  })

  output$mm30_pathway_pvals_combined_table <- renderDataTable({
    message("mm30_pathway_pvals_combined_table")

    dat <- mm30_pathways %>%
      select(-pathway) %>%
      select(pathway = pathway_html, everything()) %>%
      rownames_to_column("rank")

    dat$rank <- as.numeric(dat$rank)

    # construct data table
    float_cols <- c("all", "treatment_response")

    DT::datatable(dat, style = "bootstrap", escape = FALSE, rownames = FALSE, options = tableOpts) %>%
      formatSignif(columns = float_cols, digits = 3)
  })

  output$mm30_gene_coex_table <- renderDataTable({
    req(gene_coex())

    message("mm30_gene_coex_table")

    float_cols <- c("cor", "cor_abs")

    dat <- gene_coex()

    # construct data table
    DT::datatable(dat, style = "bootstrap", escape = FALSE, options = tableOpts) %>%
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
        print(get_feat_plot())
        dev.off()
      }
  )

  output$download_plot_tiff <- downloadHandler(
      filename = function() {
        pheno <- rv$plot_covariate
        pheno <- sub("_pval", "", pheno)

        sprintf("%s.tiff", pheno)
      },

      content = function(file) {
        tiff(file, width = 6, height = 5, res = 300, units = "in")
        print(get_feat_plot())
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

    # get covariate info for selected feature-phenotype association
    covariate_info <- covariates[[dataset_id]][[covariate]]
    assoc_method <- covariate_info$method

    dat <- feature_plot_dat()

    if (assoc_method == "survival") {
      # survival plot
      message("1) plot_survival")
      plot_survival(dat, dataset_id, covariate, rv$plot_feature,
                    input$survival_cutoff_percentile, cfg$colors)
    } else if (assoc_method == "logit") {
      # violin plot
      message("2) plot_categorical")
      plot_categorical(dat, dataset_id, covariate, rv$plot_feature, cfg$colors)
    } else if (assoc_method == "deseq") {
      # scatter plot
      message("3) plot_deseq")
      plot_deseq(dat, dataset_id, covariate, rv$plot_feature, cfg$colors)
    }
  })
  output$feature_plot <- renderPlot(featurePlot())

  get_feat_plot <- function() {
    featurePlot()
  }

  # co-expression tab
  output$gene_coex_plot <- renderPlot({
    req(mm30_expr)
    req(input$coex_gene1)
    req(input$coex_gene2)

    gene1 <- input$coex_gene1
    gene2 <- input$coex_gene2

    dat_name <- input$coex_experiment

    dat <- mm30_expr[[dat_name]] %>%
      filter(symbol %in% c(gene1, gene2))

    if (nrow(dat) == 2) {
      dat_long <- dat %>%
        column_to_rownames("symbol") %>%
        t() %>%
        as.data.frame()

      # compute pearson correlation
      coex_cor <- cor(dat_long[, 1], dat_long[, 2])

      ggplot(dat_long, aes_string(gene1, gene2)) +
          geom_point() +
          geom_smooth(method = lm) +
          ggtitle(sprintf("%s (cor = %0.2f)", dat_name, coex_cor))
    } else {
      # if both genes cannot be found in the dataset, for now, don't do anything..
      NULL
    }
  })

  # bookmarking support
  observe({
    reactiveValuesToList(input)
    session$doBookmark()
  })
  onBookmarked(updateQueryString)

  setBookmarkExclude(c(
    "mm30_gene_pvals_combined_table_rows_selected",
    "mm30_gene_pvals_combined_table_columns_selected",
    "mm30_gene_pvals_combined_table_cells_selected",
    "mm30_gene_pvals_combined_table_rows_current",
    "mm30_gene_pvals_combined_table_rows_all",
    "mm30_gene_pvals_combined_table_state",
    "mm30_gene_pvals_combined_table_search",
    "mm30_gene_pvals_combined_table_cell_clicked",
    "mm30_pathway_pvals_combined_table_rows_selected",
    "mm30_pathway_pvals_combined_table_columns_selected",
    "mm30_pathway_pvals_combined_table_cells_selected",
    "mm30_pathway_pvals_combined_table_rows_current",
    "mm30_pathway_pvals_combined_table_rows_all",
    "mm30_pathway_pvals_combined_table_state",
    "mm30_pathway_pvals_combined_table_search",
    "mm30_pathway_pvals_combined_table_cell_clicked",
    "mm30_gene_coex_table_rows_selected",
    "mm30_gene_coex_table_columns_selected",
    "mm30_gene_coex_table_cells_selected",
    "mm30_gene_coex_table_rows_current",
    "mm30_gene_coex_table_rows_all",
    "mm30_gene_coex_table_state",
    "mm30_gene_coex_table_search",
    "mm30_gene_coex_table_cell_clicked"
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
      windowTitle = "mm30",

      tabPanel(
        "Genes",
        conditionalPanel(
          condition = "typeof input.mm30_gene_pvals_combined_table_rows_selected  !== 'undefined' 
                           && input.mm30_gene_pvals_combined_table_rows_selected.length > 0",
          column(
            width = 4,
            uiOutput("gene_info")
          )
        ),
        withSpinner(dataTableOutput("mm30_gene_pvals_combined_table"))
      ),
      tabPanel(
        "Pathways",
        withSpinner(dataTableOutput("mm30_pathway_pvals_combined_table"))
      ),
      tabPanel(
        "Visualize",
        fluidRow(
          column(
            width = 4,
            selectInput("feature_type", "Feature Level:", choices = feature_types,
                        selected = "genes"),
            helpText("Feature level to visualize."),
            hr(),
            selectizeInput("feature", "Feature:", choices = NULL),
            helpText("Feature to visualize."),
            hr(),
            uiOutput("plot_covariate"),
            helpText("Phenotype / covariate to compare feature expression or SNP counts against."),
            hr(),
            selectInput("survival_cutoff_percentile", "Upper/Lower Feature Expression Cutoffs:",
                        choices = surv_expr_cutoffs, selected = 25),
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
          selectizeInput("coex_gene1", "Gene 1", choices = NULL),
          selectizeInput("coex_gene2", "Gene 2", choices = NULL),
          selectizeInput("coex_experiment", "Experiment", choices = NULL),
          hr(),
          helpText("If present, the table lists the top 100 most highly co-expressed genes with ",
                   "\"gene1\" in the selected dataset.",
                   "At present, only a subset of the top ~25% of genes with the highest variance ",
                   "are included in the analysis; In the future, this will be expanded to ",
                   "include all genes present in the dataset."),

        ),
        column(
          width = 3,
          withSpinner(plotOutput("gene_coex_plot", height = "740px"))
        ),
        column(
          width = 3,
          withSpinner(dataTableOutput("mm30_gene_coex_table"))
        )
      )
    )
  )
}

shinyApp(ui = ui, server = server, enableBookmarking = "url")
