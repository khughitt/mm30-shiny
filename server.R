#############################
# Shiny Server
#############################
server <- function(input, output, session) {
  log_info("shiny::server()")

  ############################################################ 
  # Overall survival
  ############################################################ 
  surv_os_datasets <- covariate_mdata %>% 
    filter(phenotype=='overall_survival') %>% 
    pull(dataset) %>%
    sort(TRUE)

  surv_os_gene_selected <- reactive({
    surv_os_gene_scores$Gene[input$surv_os_gene_tbl_rows_selected]
  })

  surv_os_pathway_selected <- reactive({
    surv_os_pathway_scores$Pathway[input$surv_os_pathway_tbl_rows_selected]
  })

  surv_os_gene_details <- reactive({
    req(input$surv_os_gene_tbl_rows_selected)

    gene <- surv_os_gene_selected()

    df <- bind_rows(
      surv_os_gene_pvals %>%
        filter(symbol == gene),
      surv_os_gene_effects %>%
        filter(symbol == gene),
      surv_os_gene_errors %>%
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

  surv_os_pathway_details <- reactive({
    req(input$surv_os_pathway_tbl_rows_selected)

    pathway <- surv_os_pathway_selected()

    df <- bind_rows(
      surv_os_pathway_pvals %>%
        filter(gene_set == pathway),
      surv_os_pathway_effects %>%
        filter(gene_set == pathway),
      surv_os_pathway_errors %>%
        filter(gene_set == pathway)
    )

    df %>%
      select(-gene_set) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("P-value", "Hazard Ratio", "Std. Error")) %>%
      rownames_to_column("Dataset") %>%
      mutate(Dataset=str_replace(Dataset, '_overall_survival', ''))
  })

  surv_os_gene_tcga <- reactive({
    req(input$surv_os_gene_tbl_rows_selected)
    gene <- surv_os_gene_selected()

    tcga %>%
      filter(symbol == gene) %>%
      select(-symbol, -ensgene)
  })

  output$surv_os_gene_summary_html <- renderUI({
      req(input$surv_os_gene_tbl_rows_selected)
      gene <- surv_os_gene_selected()

      all_rank <- gene_scores_all %>% 
        filter(symbol == gene) %>%
        pull(Rank)

      surv_os_rank <- surv_os_gene_scores %>% 
        filter(Gene == gene) %>%
        pull(Rank)

      surv_pfs_rank <- gene_scores_surv_pfs %>% 
        filter(Gene == gene) %>%
        pull(Rank)

      tagList(
        tags$h2(gene),
        br(),
        tags$h4("MM30 Rankings:"),
        br(),
        HTML(sprintf("All: <b>%d</b>",  all_rank)), 
        br(),
        HTML(sprintf("Overall Survival (OS): <b>%d</b>", surv_os_rank)), 
        br(),
        HTML(sprintf("Progression-free Survival (PFS): <b>%d</b>", surv_pfs_rank))
      )
  })

  output$surv_os_pathway_summary_html <- renderUI({
    tagList(tags$h2("PLACEHOLDER"))
  })

  output$surv_os_gene_details_tbl <- renderTable({
    surv_os_gene_details()
  }, digits=-3)

  output$surv_os_pathway_details_tbl <- renderTable({
    surv_os_pathway_details()
  }, digits=-3)

  output$surv_os_gene_tcga_tbl <- renderTable({
    surv_os_gene_tcga()
  }, digits=4)

  output$surv_os_gene_tbl <- renderDT({
    DT::datatable(surv_os_gene_scores, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$surv_os_pathway_tbl <- renderDT({
    DT::datatable(surv_os_pathway_scores, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$surv_os_plot <- renderPlot({
    req(input$surv_os_gene_tbl_rows_selected)

    plts <- list()

    # lapply(surv_os_datasets, function(acc) {
    for (acc in surv_os_datasets) {

      df <- indiv_mdata[[acc]] %>%
        select(1, time=os_time, event=os_censor)

      sample_ids <- df %>%
        pull(1)

      gene_name <- surv_os_gene_selected()
      gene_expr <- as.numeric(expr_dat[gene_name, sample_ids])

      df$feature <- gene_expr 

      time_units <- "Days"

      if (acc %in% c("GSE106218", "GSE19784", "GSE24080")) {
        time_units <- "Months"
      }

      if (sum(!is.na(gene_expr)) > 0) {
        plts[[acc]] <- plot_survival(df, acc, time_units, surv_expr_cutoffs, cfg$colors)
      }

    # })
    }

    arrange_ggsurvplots(plts, nrow=ceiling(length(plts) / 3), ncol=3,
                        title=sprintf("Overall survival (%s)", gene_name))
  })

  ##########
  #
  # Cell lines (HMCL)
  #
  ###########
  output$hmcl_expr_tbl <- renderDT({
    DT::datatable(expr_hmcl, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=hmcl_cells, digits=3)
  })

  hmcl_gene_selected <- reactive({
    expr_hmcl$symbol[input$hmcl_expr_tbl_rows_selected]
  })

  output$hmcl_hist_plot <- renderPlot({
    req(input$hmcl_expr_tbl_rows_selected)

    gene <- hmcl_gene_selected()

    x <- expr_hmcl %>%
      filter(symbol == gene) %>%
      select(-symbol) %>%
      as.numeric()

    df <- data.frame(x=x)

    ggplot(df, aes(x)) +
      geom_histogram() +
      ggtitle(sprintf("HMCL Expression (%s)", gene)) +
      xlab("Expression (TPM)") +
      ylab("Count")
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
