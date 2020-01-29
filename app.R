#
# Feature weights shiny UI
#
library(DT)
library(feather)
library(ggdark)
library(gridExtra)
library(plotly)
library(shinythemes)
library(tidyverse)
library(yaml)

options(stringsAsFactors = FALSE)
options(digits = 3)

# load feature/response data
# rna <- read_feather('/data/nih/pharmacogx/1.0/gdsc1000/features/rna.feather')
# auc <- read_feather('/data/nih/pharmacogx/1.0/gdsc1000/phenotypes/auc_recomputed.feather')

# load MM25 combined scores
mm25_genes <- read_feather('/data/nih/fassoc/1.0/summary/gene_combined_scores.feather')
mm25_pathways <- read_feather('/data/nih/fassoc/1.0/summary/pathway_combined_scores.feather')

mm25_gene_pvals <- read_feather('/data/nih/fassoc/1.0/summary/gene_scores.feather')

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

server <- function(input, output, session) {
  # reactives
  gene_pvals <- reactive({
    req(input$select_gene)

    dat <- tibble(
      phenotype = colnames(mm25_gene_pvals)[-1],
      pval = as.numeric(mm25_gene_pvals[mm25_gene_pvals$symbol == input$select_gene, -1])
    )

    dat <- dat %>%
      arrange(pval)
  })

  output$gene_plots <- renderPlot({
    plts <- list()

    for (pheno in head(gene_pvals()$phenotype, 4)) {
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

      cat('Foo:\n')
      cat(mdat)
      cat('\n')
      cat(covariate)
      cat('\n')

      # config
      cfg <- cfgs[[dataset]]$phenotypes[[mdat]]$associations[[covariate]]

      cat('Okay...\n')

      if (cfg$method == 'survival') {
        key <- cfg$params$time
      } else {
        key <- cfg$params$field
      }

      # feature data
      feature <- gene_data[[dataset]] %>%
        filter(symbol == input$select_gene) %>%
        select(-symbol) %>%
        as.numeric()

      # pheno data
      response <- pull(pheno_data[[dataset]][[mdat]], key)

      dat <- data.frame(feature, response)

      if (cfg$method == 'logit') {
        response <- factor(response)

        plt <- ggplot(dat, aes(x = response, y = feature)) +
          geom_boxplot() +
          ggtitle(sprintf("%s: %s vs. %s", dataset, input$select_gene, covariate)) +
          xlab(covariate) +
          ylab(sprintf("%s expression", input$select_gene)) +
          dark_theme_gray()
      } else {
        plt <- ggplot(dat, aes(x = feature, y = response)) +
          geom_point(alpha = 0.75) +
          ggtitle(sprintf("%s: %s vs. %s", dataset, covariate, input$select_gene)) +
          xlab(sprintf("%s expression", input$select_gene)) +
          ylab(covariate) +
          dark_theme_gray()
      }

      plts[[pheno]] <- plt
    }

    print(grid.arrange(grobs = plts, ncol = 2))
  })

  # table output
  output$mm25_genes <- renderDataTable({
    DT::datatable(mm25_genes, style = 'bootstrap', options = list(pageLength = 15)) %>%
      formatRound(columns = float_cols, digits = 5)
  })

  output$mm25_pathways <- renderDataTable({
    DT::datatable(mm25_pathways, style = 'bootstrap', options = list(pageLength = 15)) %>%
      formatRound(columns = float_cols, digits = 5)
  })

  # form output
  # output$select_gene <- renderUI({
  #   selectizeInput('select_gene', "Gene:", choices = as.character(mm25_genes$symbol)[1:10],
  #               selected = as.character(mm25_genes$symbol)[1])
  # })
  updateSelectizeInput(session, 'select_gene', 
                       choices = as.character(mm25_genes$symbol), server = TRUE)

  # plot output
  # output$weights_dist <- renderPlotly({
  #   scores <- weights_table() %>%
  #     pull(score)
  #
  #   dat <- data.frame(score = scores)
  #
  #   ggplot(dat, aes(x = score)) +
  #     geom_density(alpha = 0.85, fill = '#9dff00') +
  #     ggtitle(sprintf("Distribution of %s", input$score_field)) +
  #     dark_theme_gray() +
  #     theme(
  #       plot.background = element_rect(fill = "#222222"),
  #       panel.grid.major = element_line(color = "#555555", size = 0.2),
  #       panel.grid.minor = element_line(color = "#555555", size = 0.2)
  #     )
  # })
}

ui <- fluidPage(
  titlePanel("MM25 (Jan 29, 2020)", windowTitle = "MM25 (Jan 29, 2020)"),
  tags$head(includeCSS("resources/styles.css")),
  theme = shinytheme("darkly"),
  tabsetPanel(
    tabPanel(
      "Genes",
      dataTableOutput('mm25_genes')
    ),
    tabPanel(
      "Pathways",
      dataTableOutput('mm25_pathways')
    ),
    tabPanel(
      "Datasets",
      column(
        width = 3,
        selectizeInput('select_gene', "Gene:", choices = NULL)
      ),
      column(
        width = 9,
        plotOutput('gene_plots', height = '920px')
      )
    )
  )
)

shinyApp(ui = ui, server = server)
