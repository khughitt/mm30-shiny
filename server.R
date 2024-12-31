#############################
# Shiny Server
#############################
server <- function(input, output, session) {
  log_info("shiny::server()")

  ############################################################
  # Overall survival
  ############################################################

  #---------------------------------------
  # Overall survival > genes
  #---------------------------------------
  surv_os_datasets <- covariate_mdata %>%
    filter(phenotype=='overall_survival') %>%
    pull(dataset) %>%
    sort(TRUE)

  surv_os_gene_selected <- reactive({
    surv_os_gene_scores$symbol[input$surv_os_gene_scores_tbl_rows_selected]
  })

  surv_os_gene_details <- reactive({
    req(input$surv_os_gene_scores_tbl_rows_selected)

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

  tcga_associations_tbl <- reactive({
    req(!is.null(selectedGene()))
    gene <- selectedGene()

    tcga %>%
        filter(symbol == gene) %>%
        select(-symbol, -ensgene)
  })

  output$surv_os_gene_details_tbl <- renderTable({
    surv_os_gene_details()
  }, digits=-3)

  output$surv_os_gene_tcga_tbl <- output$stage_gene_tcga_tbl <- output$treatment_gene_tcga_tbl <- renderTable({
    # rv$tcga_tbl
    tcga_associations_tbl()
  }, digits=4)

  output$surv_os_gene_scores_tbl <- renderDT({
    df <- surv_os_gene_scores %>%
      select(Gene=symbol, Rank, `P-value\n(metap)`=sumz_wt_pval,
             `P-value\n(metafor)`=metafor_pval, Description=description,
             CHR=chr_region, `Cell Cycle`=cell_cycle_phase, `DGIdb\ncategories`=dgidb_categories)

    DT::datatable(df, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$surv_os_gene_plot <- renderPlot({
    req(input$surv_os_gene_scores_tbl_rows_selected)

    plts <- list()

    for (acc in surv_os_datasets) {

      df <- indiv_mdata[[acc]] %>%
        select(1, time=os_time, event=os_censor)

      sample_ids <- df %>%
        pull(1)

      gene_name <- surv_os_gene_selected()
      feat_expr <- as.numeric(gene_expr[gene_name, sample_ids])

      df$feature <- feat_expr

      time_units <- "Days"

      if (acc %in% c("GSE106218", "GSE19784", "GSE24080")) {
        time_units <- "Months"
      }

      if (sum(!is.na(feat_expr)) > 0) {
        plts[[acc]] <- plot_survival(df, acc, time_units, surv_expr_cutoffs, cfg$colors)
      }
    }

    arrange_ggsurvplots(plts, nrow=ceiling(length(plts) / 3), ncol=3,
                        title=sprintf("Overall survival (%s)", gene_name))
  })

  #---------------------------------------
  # Overall survival > gene sets
  #---------------------------------------
  surv_os_gene_set_selected <- reactive({
    surv_os_gene_set_scores$gene_set[input$surv_os_gene_set_scores_tbl_rows_selected]
  })

  surv_os_gene_set_details <- reactive({
    req(input$surv_os_gene_set_scores_tbl_rows_selected)

    selected_gene_set <- surv_os_gene_set_selected()

    df <- bind_rows(
      surv_os_gene_set_pvals %>%
        filter(gene_set == selected_gene_set),
      surv_os_gene_set_effects %>%
        filter(gene_set == selected_gene_set),
      surv_os_gene_set_errors %>%
        filter(gene_set == selected_gene_set)
    )

    df %>%
      select(-gene_set) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("P-value", "Hazard Ratio", "Std. Error")) %>%
      rownames_to_column("Dataset") %>%
      mutate(Dataset=str_replace(Dataset, '_overall_survival', ''))
  })

  output$surv_os_gene_set_details_tbl <- renderTable({
    surv_os_gene_set_details()
  }, digits=-3)

  output$surv_os_gene_set_gene_tbl <- renderTable({
    req(input$surv_os_gene_set_scores_tbl_rows_selected)

    #gene_set <- surv_os_gene_set_scores$gene_set[input$surv_os_gene_set_scores_tbl_rows_selected]
    selected <- surv_os_gene_set_selected()

    genes <- gene_set_mdata %>%
      filter(gene_set == selected) %>%
      pull(genes) %>%
      unlist()

    surv_os_gene_scores %>%
      filter(symbol %in% genes) %>%
      select(Gene=symbol, `Rank (OS)`=Rank, `P-value\n(metap)`=sumz_wt_pval, `P-value (metafor)`=metafor_pval,
             Description=description, CHR=chr_region, `Cell Cycle`=cell_cycle_phase,
             `DGIdb\ncategories`=dgidb_categories)

  }, digits=4)

  output$surv_os_gene_set_scores_tbl <- renderDT({
    DT::datatable(surv_os_gene_set_scores, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$surv_os_gene_set_plot <- renderPlot({
    req(input$surv_os_gene_set_scores_tbl_rows_selected)

    plts <- list()

    for (acc in surv_os_datasets) {

      df <- indiv_mdata[[acc]] %>%
        select(1, time=os_time, event=os_censor)

      sample_ids <- df %>%
        pull(1)

      selected_gene_set <- surv_os_gene_set_selected()

      feat_expr <- as.numeric(gene_set_expr[selected_gene_set, sample_ids])

      df$feature <- feat_expr

      time_units <- "Days"

      if (acc %in% c("GSE106218", "GSE19784", "GSE24080")) {
        time_units <- "Months"
      }

      if (sum(!is.na(feat_expr)) > 0) {
        plts[[acc]] <- plot_survival(df, acc, time_units, surv_expr_cutoffs, cfg$colors)
      }
    }

    arrange_ggsurvplots(plts, nrow=ceiling(length(plts) / 3), ncol=3,
                        title=sprintf("Overall survival (%s)", selected_gene_set))
  })

  ############################################################
  # Disease Stage
  ############################################################
  stage_datasets <- covariate_mdata %>%
    filter(phenotype=='disease_stage') %>%
    pull(dataset) %>%
    sort(TRUE)

  stage_gene_selected <- reactive({
    stage_gene_scores$symbol[input$stage_gene_scores_tbl_rows_selected]
  })

  stage_gene_set_selected <- reactive({
    stage_gene_set_scores$gene_set[input$stage_gene_set_scores_tbl_rows_selected]
  })

  stage_gene_details <- reactive({
    req(input$stage_gene_scores_tbl_rows_selected)

    gene <- stage_gene_selected()

    df <- bind_rows(
      stage_gene_pvals %>%
        filter(symbol == gene),
      stage_gene_effects %>%
        filter(symbol == gene),
      stage_gene_errors %>%
        filter(symbol == gene)
    )

    df %>%
      select(-symbol) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("P-value", "Odds Ratio", "Std. Error")) %>%
      rownames_to_column("Dataset") %>%
      mutate(Dataset=str_replace(Dataset, '_disease_stage', ''))
  })

  stage_gene_set_details <- reactive({
    req(input$stage_gene_set_scores_tbl_rows_selected)

    selected_gene_set <- stage_gene_set_selected()

    df <- bind_rows(
      stage_gene_set_pvals %>%
        filter(gene_set == selected_gene_set),
      stage_gene_set_effects %>%
        filter(gene_set == selected_gene_set),
      stage_gene_set_errors %>%
        filter(gene_set == selected_gene_set)
    )

    df %>%
      select(-gene_set) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("P-value", "Odds Ratio", "Std. Error")) %>%
      rownames_to_column("Dataset") %>%
      mutate(Dataset=str_replace(Dataset, '_disease_stage', ''))
  })

  output$stage_gene_summary_html <- output$surv_os_gene_summary_html <- output$surv_pfs_gene_summary_html <- renderUI({
      req(!is.null(selectedGene()))
      gene <- selectedGene()

      all_rank <- gene_scores_all %>%
        filter(symbol == gene) %>%
        pull(Rank)

      surv_os_rank <- surv_os_gene_scores %>%
        filter(symbol == gene) %>%
        pull(Rank)

      surv_pfs_rank <- surv_pfs_gene_scores %>%
        filter(symbol == gene) %>%
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

  output$stage_gene_set_summary_html <- output$surv_os_gene_set_summary_html  <- output$surv_pfs_gene_set_summary_html <- renderUI({
      req(!is.null(selectedGeneSet()))
      selected_gene_set <- selectedGeneSet()

      all_rank <- gene_set_scores_all %>%
        filter(gene_set == selected_gene_set) %>%
        pull(Rank)

      surv_os_rank <- surv_os_gene_set_scores %>%
        filter(gene_set == selected_gene_set) %>%
        pull(Rank)

      surv_pfs_rank <- surv_pfs_gene_set_scores %>%
        filter(gene_set == selected_gene_set) %>%
        pull(Rank)

      gene_set_info <- gene_set_mdata %>%
        filter(gene_set == selected_gene_set)

      collection <- gene_set_info %>%
        pull(collection)

      url <- gene_set_info %>%
        pull(url)

      num_genes <- gene_set_info %>%
        pull(collection_size)

      tagList(
        tags$h2(sprintf("%s / %s (n=%d)", collection, selected_gene_set, num_genes)),
        br(),
        tags$a(url, target="_blank", href=url),
        tags$h4("MM30 Rankings:"),
        br(),
        HTML(sprintf("All: <b>%d</b>",  all_rank)),
        br(),
        HTML(sprintf("Overall Survival (OS): <b>%d</b>", surv_os_rank)),
        br(),
        HTML(sprintf("Progression-free Survival (PFS): <b>%d</b>", surv_pfs_rank))
      )
  })


  output$stage_gene_details_tbl <- renderTable({
    stage_gene_details()
  }, digits=-3)

  output$stage_gene_set_details_tbl <- renderTable({
    stage_gene_set_details()
  }, digits=-3)

  output$stage_gene_set_gene_tbl <- renderTable({
    req(input$stage_gene_set_scores_tbl_rows_selected)

    selected_gene_set <- stage_gene_set_scores$gene_set[input$stage_gene_set_scores_tbl_rows_selected]

    genes <- gene_set_mdata %>%
      filter(gene_set == selected_gene_set) %>%
      pull(genes) %>%
      unlist()

    stage_gene_scores %>%
      filter(symbol %in% genes) %>%
      select(Gene=symbol, `Rank (OS)`=Rank, `P-value\n(metap)`=sumz_wt_pval, `P-value (metafor)`=metafor_pval,
             Description=description, CHR=chr_region, `Cell Cycle`=cell_cycle_phase,
             `DGIdb\ncategories`=dgidb_categories)

  }, digits=4)

  output$stage_gene_scores_tbl <- renderDT({
    df <- stage_gene_scores %>%
      select(Gene=symbol, Rank, `P-value\n(metap)`=sumz_wt_pval,
             `P-value\n(metafor)`=metafor_pval, Description=description,
             CHR=chr_region, `Cell Cycle`=cell_cycle_phase, `DGIdb\ncategories`=dgidb_categories)

    DT::datatable(df, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$stage_gene_set_scores_tbl <- renderDT({
    df <- stage_gene_set_scores %>%
      select(`Gene set`=gene_set, everything())

    DT::datatable(df, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  output$stage_gene_set_plot <- renderPlot({
    req(input$stage_gene_set_scores_tbl_rows_selected)

    gene_set_name <- stage_gene_set_selected()

    df <- stage_gene_set_scaled_expr %>%
      filter(gene_set == gene_set_name) %>%
      select(-gene_set)

    ggplot(df, aes(x=stage, y=expr, fill=stage)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=cfg$colors) +
      facet_wrap(~dataset, ncol=3, scales='free')
  })

  output$stage_gene_plot <- renderPlot({
    req(input$stage_gene_scores_tbl_rows_selected)

    gene_name <- stage_gene_selected()

    df <- stage_gene_scaled_expr %>%
      filter(symbol == gene_name) %>%
      select(-symbol)

    ggplot(df, aes(x=stage, y=expr, fill=stage)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=cfg$colors) +
      facet_wrap(~dataset, ncol=3, scales='free')
  })

  ############################################################
  # Treatment response
  ############################################################

  #---------------------------------------
  # Treatment response > genes
  #---------------------------------------
  treatment_datasets <- covariate_mdata %>%
    filter(phenotype=='treatment_response') %>%
    pull(dataset) %>%
    sort(TRUE)

  treatment_gene_selected <- reactive({
    treatment_gene_scores$symbol[input$treatment_gene_scores_tbl_rows_selected]
  })

  output$treatment_gene_scores_tbl <- renderDT({
    df <- treatment_gene_scores %>%
      select(Gene=symbol, Rank, `P-value\n(metap)`=sumz_wt_pval,
             `P-value\n(metafor)`=metafor_pval, Description=description,
             CHR=chr_region, `Cell Cycle`=cell_cycle_phase, `DGIdb\ncategories`=dgidb_categories)

    DT::datatable(df, style="bootstrap", escape=FALSE,  rownames=FALSE,
                  selection=selectOpts, options=tableOpts) %>%
      formatSignif(columns=c("P-value\n(metap)", "P-value\n(metafor)"), digits=3)
  })

  treatment_gene_details <- reactive({
    req(input$treatment_gene_scores_tbl_rows_selected)

    gene <- treatment_gene_selected()

    df <- bind_rows(
      treatment_gene_pvals %>%
        filter(symbol == gene),
      treatment_gene_effects %>%
        filter(symbol == gene),
      treatment_gene_errors %>%
        filter(symbol == gene)
    )

    df %>%
      select(-symbol) %>%
      t() %>%
      as.data.frame() %>%
      setNames(c("P-value", "Effect", "Std. Error")) %>%
      rownames_to_column("Dataset") %>%
      mutate(Dataset=str_replace(Dataset, '_treatment_response', ''))
  })

  output$treatment_gene_details_tbl <- renderTable({
    treatment_gene_details()
  }, digits=-3)

  output$treatment_gene_plot <- renderPlotly({
    req(input$treatment_gene_scores_tbl_rows_selected)

    plts <- list()

    # treatment plot titles
    plot_titles <- list(
      GSE9782 = "Treatment Response (VTD)",
      GSE68871 = "Treatment Response (VTD)",
      GSE39754 = "Treatment Response (VAD + ACST)"
    )

    gene_name <- treatment_gene_selected()

    for (acc in treatment_datasets) {

      df <- indiv_mdata[[acc]] %>%
        select(1, response=treatment_response)

      sample_ids <- df %>%
        pull(1)

      feat_expr <- as.numeric(gene_expr[gene_name, sample_ids])

      df$feature <- feat_expr

      if (sum(!is.na(feat_expr)) > 0) {
        log_info(acc)
        plts[[acc]] <- plot_categorical(df, acc, plot_titles[[acc]], gene_name, cfg$colors)
      }
    }

    # MMRF plots
    df <- indiv_mdata[["MMRF"]] %>%
      select(1, response=response_bor_len_dex)

    sample_ids <- df %>%
      pull(1)

    df$feature <- as.numeric(gene_expr[gene_name, sample_ids])

    if (sum(!is.na(df$feature)) > 0) {
      plts[["mmrf1"]] <- plot_categorical(df, "MMRF", "Treatment Response (Bor-Len-Dex)", gene_name, cfg$colors)
    }

    df <- indiv_mdata[["MMRF"]] %>%
      select(1, response=response_bor_cyc_dex)

    sample_ids <- df %>%
      pull(1)

    df$feature <- as.numeric(gene_expr[gene_name, sample_ids])

    if (sum(!is.na(df$feature)) > 0) {
      plts[["mmrf2"]] <- plot_categorical(df, "MMRF", "Treatment Response (Bor-Cyc-Dex)", gene_name, cfg$colors)
    }

    # "plotly::" namespace necessary or else Hmisc::plotly may be used...
    do.call(plotly::subplot, c(plts, list(nrows=2)))
  })

  ##########
  #
  # Co-expression
  #
  ###########
  gene_coex_df <- reactive({
    feat1 <- input$coex_gene1
    feat2 <- input$coex_gene2

    df <- gene_expr[c(feat1, feat2), ] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      left_join(sample_mdata, by="sample_id") %>%
      select(-platform_id, -platform_type)


    df$disease_stage[is.na(df$disease_stage)] <- "Unknown"
    df <- df[complete.cases(df), ]

    # compute linear model r^2 for each dataset and adjust experiment field so that
    # it is included in the plot titles
    for (expt_id in unique(df$experiment)) {
      expt_df <- df %>%
        filter(experiment == expt_id)

      fit <- summary(lm(expr(!!sym(feat1) ~ !!sym(feat2)), data=expt_df))

      r2 <- fit$r.squared
      pval <- fit$coefficients[2, 'Pr(>|t|)']

      expt_label <- sprintf("%s (R^2 %0.2f)", expt_id, r2, pval)
      df$experiment[df$experiment == expt_id] <- expt_label
    }

    df
  })

  output$gene_coex_plot <- renderPlot({
    feat1 <- input$coex_gene1
    feat2 <- input$coex_gene2

    df <- gene_coex_df()

    ggplot(df, aes(x=.data[[feat1]], y=.data[[feat2]])) +
      geom_point(aes(color=disease_stage)) +
      geom_smooth(method='lm') +
      ggtitle(sprintf("Gene expression: %s vs. %s", feat2, feat1)) +
      facet_wrap(~experiment, scales='free') +
      xlab(feat1) +
      ylab(feat2)
  })

  gene_set_coex_df <- reactive({
    feat1 <- input$coex_gene_set1
    feat2 <- input$coex_gene_set2

    df <- gene_set_expr[c(feat1, feat2), ] %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      left_join(sample_mdata, by="sample_id") %>%
      select(-platform_id, -platform_type)


    df$disease_stage[is.na(df$disease_stage)] <- "Unknown"
    df <- df[complete.cases(df), ]

    # compute linear model r^2 for each dataset and adjust experiment field so that
    # it is included in the plot titles
    for (expt_id in unique(df$experiment)) {
      expt_df <- df %>%
        filter(experiment == expt_id)

      fit <- summary(lm(expr(!!sym(feat1) ~ !!sym(feat2)), data=expt_df))

      r2 <- fit$r.squared
      pval <- fit$coefficients[2, 'Pr(>|t|)']

      expt_label <- sprintf("%s (R^2 %0.2f)", expt_id, r2, pval)
      df$experiment[df$experiment == expt_id] <- expt_label
    }

    df
  })

  output$gene_set_coex_plot <- renderPlot({
    feat1 <- input$coex_gene_set1
    feat2 <- input$coex_gene_set2

    df <- gene_set_coex_df()

    ggplot(df, aes(x=.data[[feat1]], y=.data[[feat2]])) +
      geom_point(aes(color=disease_stage)) +
      geom_smooth(method='lm') +
      ggtitle(sprintf("Gene set expression: %s vs. %s", feat2, feat1)) +
      facet_wrap(~experiment, scales='free') +
      xlab(feat1) +
      ylab(feat2)
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

  # common variables to store selected gene / gene set
  selectedGene <- reactiveVal()
  selectedGeneSet <- reactiveVal()

  observeEvent(input$surv_os_gene_scores_tbl_rows_selected, {
    selectedGene(surv_os_gene_selected())
  })
  observeEvent(input$stage_gene_scores_tbl_rows_selected, {
    selectedGene(stage_gene_selected())
  })

  observeEvent(input$surv_os_gene_set_scores_tbl_rows_selected, {
    selectedGeneSet(surv_os_gene_set_selected())
  })
  observeEvent(input$stage_gene_set_scores_tbl_rows_selected, {
    selectedGeneSet(stage_gene_set_selected())
  })

  ###################################
  #
  # Datasets
  #
  ###################################
  output$datasets_tbl <- renderDT({
    df <- dataset_mdata %>%
      select(dataset, title, num_samples)

    DT::datatable(df, style="bootstrap", escape=FALSE,  rownames=FALSE, selection=selectOpts,
                  options=tableOpts)
  })

  output$dataset_summary <- renderUI({
    acc <- dataset_mdata$dataset[input$datasets_tbl_rows_selected]

    mdat <- dataset_mdata %>%
      filter(dataset==acc)  %>%
      as.list()

    tagList(
      tags$h2(mdat$title),
      br(),
      tags$a(acc, target="_blank", href=mdat$urls),
      br(),
      HTML(mdat$name),
      br(),
      HTML(sprintf("Submitted: %s", mdat$submission_date)),
      br(),
      HTML(sprintf("# samples: <b>%d</b>",  mdat$num_samples)),
      br(),
      HTML(sprintf("Covariates: <b>%s</b>",  mdat$covariates)),
      br(),
      tags$h4("Abstract:"),
      HTML(mdat$abstract),
      br(),
      tags$h4("Overall design:"),
      HTML(mdat$overall_design),
      br()
    )
  })

  # bookmarking support
  # observe({
  #   reactiveValuesToList(input)
  #   session$doBookmark()
  # })
  # onBookmarked(updateQueryString)
  #
  # setBookmarkExclude(c(
  #   "disease_stage_tbl_rows_selected",
  #   "disease_stage_tbl_columns_selected",
  #   "disease_stage_tbl_cells_selected",
  #   "disease_stage_tbl_rows_current",
  #   "disease_stage_tbl_rows_all",
  #   "disease_stage_tbl_state",
  #   "disease_stage_tbl_search",
  #   "disease_stage_tbl_cell_clicked",
  #   "treatment_response_tbl_rows_selected",
  #   "treatment_response_tbl_columns_selected",
  #   "treatment_response_tbl_cells_selected",
  #   "treatment_response_tbl_rows_current",
  #   "treatment_response_tbl_rows_all",
  #   "treatment_response_tbl_state",
  #   "treatment_response_tbl_search",
  #   "treatment_response_tbl_cell_clicked"
  # ))
}
