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
library(NMF)
library(plotly)
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(survminer)
library(svglite)
library(tidyverse)
library(viridis)
library(yaml)

source("R/plotting.R")

#############################
# Setup
#############################
options(spinner.color="#00bc8c")
options(scipen=2, digits=3)

#log_threshold(DEBUG)
log_info("Initializing MM30 shiny app..")

cfg <- read_yaml("config/config-v7.3.yml")
set.seed(cfg$random_seed)

enableBookmarking("url")

tableOpts  <- list(pageLength=20, scrollX=TRUE)
selectOpts <- list(target="row", mode="single", selected=1)

# ordered vector of disease stages
# disease_stages <- c("Healthy", "MGUS", "SMM", "MM", "RRMM")

# options for survival plot upper/lower expression cutoffs
# surv_expr_cutoffs <- cfg$surv_cutoff_opts
# names(surv_expr_cutoffs) <- paste0(surv_expr_cutoffs, " %")
surv_expr_cutoffs <- c(25, 75)

#############################
# Load data
#############################
log_info("Loading data..")

results_dir <- file.path(cfg$data_dir, "results")
fassoc_dir  <- file.path(cfg$data_dir, "fassoc")
mmrf_dir    <- file.path(cfg$data_dir, "mmrf", "IA22")
geo_dir     <- file.path(cfg$data_dir, "geo")

#--------------------------
# 1. Metadata
#--------------------------
gene_mdata <- read_feather(file.path(results_dir, "metadata/genes.feather")) %>%
  select(-chr_subband)

gene_set_mdata <- read_feather(file.path(results_dir, "metadata/gene_sets.feather")) %>%
  mutate(collection_size=lengths(genes))

sample_mdata <- read_feather(file.path(results_dir, "metadata/samples.feather")) %>%
  rename(dataset=experiment)

sample_counts <- sample_mdata %>%
  group_by(dataset) %>%
  summarize(num_samples=n())

covariate_mdata <- read_feather(file.path(results_dir, "metadata/covariates.feather")) %>%
  rename(dataset=accession)

dataset_covariates <- covariate_mdata %>%
  select(dataset, category) %>%
  distinct() %>%
  group_by(dataset) %>%
  summarise(covariates=paste(category, collapse=", "))

dataset_mdata <- read_feather(file.path(results_dir, "metadata/datasets.feather")) %>%
  rename(dataset=accession, title=geo_title) %>%
  select(-num_samples) %>%
  left_join(sample_counts, by='dataset') %>%
  left_join(dataset_covariates, by='dataset') %>%
  arrange(-num_samples)

dataset_mdata$title[dataset_mdata$dataset == "MMRF"] <- "MMRF CoMMpaSS Trial"

# create a list of sample metadata dataframes for each dataset
indiv_mdata <- list()

# create mmrf sample metadata table
infile <- file.path(mmrf_dir, "clinical_flat_files", "MMRF_CoMMpass_IA22_combined_metadata.feather")

mmrf_mdata <- read_feather(infile) %>%
      select(public_id, iss_stage, ecog, mm_status, fresp, frespcd,
              trtshnm, pfs_time, pfs_censor, os_time, os_censor)

mmrf_mdata$frespcd <- ordered(mmrf_mdata$frespcd, c("sCR", "CR", "VGPR", "PR"))
mmrf_mdata$response_bor_len_dex <- mmrf_mdata$frespcd
mmrf_mdata$response_bor_cyc_dex <- mmrf_mdata$frespcd
mmrf_mdata$response_bor_len_dex[mmrf_mdata$trtshnm != "Bor-Len-Dex"] <- NA
mmrf_mdata$response_bor_cyc_dex[mmrf_mdata$trtshnm != "Bor-Cyc-Dex"] <- NA

indiv_mdata[["MMRF"]] <- mmrf_mdata

# add geo sample metadata tables
for (dir_ in Sys.glob(file.path(geo_dir, "final", "GSE*"))) {
  acc <- basename(dir_)
  indiv_mdata[[acc]] <- read_feather(file.path(dir_, 'column-metadata.feather')) %>%
    ungroup()
}

#--------------------------------
# 2. "All"
#--------------------------------
gene_scores_all <- read_feather(file.path(results_dir, "scores/metap/all/gene.feather")) %>%
  left_join(gene_mdata, by='symbol')

gene_scores_all$Rank <- 1:nrow(gene_scores_all)

gene_set_scores_all <- read_feather(file.path(results_dir, "scores/metap/all/gene_set.feather")) %>%
  left_join(gene_set_mdata, by='gene_set')

gene_set_scores_all$Rank <- 1:nrow(gene_set_scores_all)

#--------------------------------
# 3. Overall survival
#--------------------------------
surv_os_gene_scores <- read_feather(file.path(results_dir, "scores/combined/survival_os/gene.feather")) %>%
    left_join(gene_mdata, by='symbol') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$surv_os) %>%
    select(symbol, sumz_wt_pval, metafor_pval, description,
           chr_region, cell_cycle_phase, dgidb_categories) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    select(Rank, everything())

surv_os_gene_pvals   <- read_feather(file.path(results_dir, "associations/survival_os/gene/pvals.feather"))
surv_os_gene_errors  <- read_feather(file.path(results_dir, "associations/survival_os/gene/errors.feather"))
surv_os_gene_effects <- read_feather(file.path(results_dir, "associations/survival_os/gene/effects.feather"))

surv_os_gene_set_scores <- read_feather(file.path(results_dir, "scores/combined/survival_os/gene_set.feather")) %>%
    left_join(gene_set_mdata, by='gene_set') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$surv_os) %>%
    select(gene_set, sumz_wt_pval,
           `P-value\n(metafor)`=metafor_pval, Collection=collection,
           `# Genes`=collection_size) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

surv_os_gene_set_pvals   <- read_feather(file.path(results_dir, "associations/survival_os/gene_set/pvals.feather"))
surv_os_gene_set_errors  <- read_feather(file.path(results_dir, "associations/survival_os/gene_set/errors.feather"))
surv_os_gene_set_effects <- read_feather(file.path(results_dir, "associations/survival_os/gene_set/effects.feather"))

#--------------------------------
# 4. Progression free survival
#--------------------------------
surv_pfs_gene_scores <- read_feather(file.path(results_dir, "scores/combined/survival_pfs/gene.feather")) %>%
    left_join(gene_mdata, by='symbol') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$surv_pfs) %>%
    select(symbol, sumz_wt_pval, `P-value\n(metafor)`=metafor_pval, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase, `DGIdb\ncategories`=dgidb_categories) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

surv_pfs_gene_set_scores <- read_feather(file.path(results_dir, "scores/combined/survival_pfs/gene_set.feather")) %>%
    left_join(gene_set_mdata, by='gene_set') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$surv_pfs) %>%
    select(gene_set, sumz_wt_pval,
           `P-value\n(metafor)`=metafor_pval, Collection=collection,
           `# Genes`=collection_size) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

surv_pfs_gene_set_pvals   <- read_feather(file.path(results_dir, "associations/survival_pfs/gene_set/pvals.feather"))
surv_pfs_gene_set_errors  <- read_feather(file.path(results_dir, "associations/survival_pfs/gene_set/errors.feather"))
surv_pfs_gene_set_effects <- read_feather(file.path(results_dir, "associations/survival_pfs/gene_set/effects.feather"))

#--------------------------------
# x. Disease stage
#--------------------------------
gene_stage_scores <- read_feather(file.path(results_dir, "scores/combined/disease_stage/gene.feather")) %>%
    left_join(gene_mdata, by='symbol') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$disease_stage) %>%
    select(symbol, sumz_wt_pval, metafor_pval, description,
           chr_region, cell_cycle_phase, dgidb_categories) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    select(Rank, everything())

gene_stage_pvals   <- read_feather(file.path(results_dir, "associations/disease_stage/gene/pvals.feather"))
gene_stage_errors  <- read_feather(file.path(results_dir, "associations/disease_stage/gene/errors.feather"))
gene_stage_effects <- read_feather(file.path(results_dir, "associations/disease_stage/gene/effects.feather"))

gene_stage_scaled_expr <- read_feather(file.path(results_dir, "disease_stage/scaled/gene", "combined.feather")) %>%
  filter(!stage %in% c('early', 'late', 'pre_relapsed'))

ind <- match(gene_stage_scaled_expr$dataset, dataset_mdata$dataset)
gene_stage_scaled_expr$dataset_name <- dataset_mdata$name[ind]

gene_set_stage_scores <- read_feather(file.path(results_dir, "scores/combined/disease_stage/gene_set.feather")) %>%
    left_join(gene_set_mdata, by='gene_set') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$disease_stage) %>%
    select(gene_set, sumz_wt_pval,
           `P-value\n(metafor)`=metafor_pval, Collection=collection,
           `# Genes`=collection_size) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

gene_set_stage_pvals   <- read_feather(file.path(results_dir, "associations/disease_stage/gene_set/pvals.feather"))
gene_set_stage_errors  <- read_feather(file.path(results_dir, "associations/disease_stage/gene_set/errors.feather"))
gene_set_stage_effects <- read_feather(file.path(results_dir, "associations/disease_stage/gene_set/effects.feather"))

gene_set_stage_scaled_expr <- read_feather(file.path(results_dir, "disease_stage/scaled/gene_set", "combined.feather")) %>%
  filter(!stage %in% c('early', 'late', 'pre_relapsed'))

ind <- match(gene_set_stage_scaled_expr$dataset, dataset_mdata$dataset)
gene_set_stage_scaled_expr$dataset_name <- dataset_mdata$name[ind]

transitions_dir <- file.path(results_dir, "disease_stage/transitions")
transition_files <- Sys.glob(file.path(transitions_dir, "gene/*.feather"))
transitions <- sub(".feather", "", basename(transition_files))

gene_transitions <- lapply(transition_files, read_feather)
names(gene_transitions) <- transitions

#--------------------------------
# x. Treatment
#--------------------------------
treatment_gene_scores <- read_feather(file.path(results_dir, "scores/combined/treatment_response/gene.feather")) %>%
    left_join(gene_mdata, by='symbol') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$treatment_response) %>%
    select(symbol, sumz_wt_pval, metafor_pval, description,
           chr_region, cell_cycle_phase, dgidb_categories) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    select(Rank, everything())

treatment_gene_pvals   <- read_feather(file.path(results_dir, "associations/treatment_response/gene/pvals.feather"))
treatment_gene_errors  <- read_feather(file.path(results_dir, "associations/treatment_response/gene/errors.feather"))
treatment_gene_effects <- read_feather(file.path(results_dir, "associations/treatment_response/gene/effects.feather"))

treatment_gene_set_scores <- read_feather(file.path(results_dir, "scores/combined/treatment_response/gene_set.feather")) %>%
    left_join(gene_set_mdata, by='gene_set') %>%
    filter(num_missing <= cfg$scores_max_missing_datasets$treatment_response) %>%
    select(gene_set, sumz_wt_pval,
           `P-value\n(metafor)`=metafor_pval, Collection=collection,
           `# Genes`=collection_size) %>%
    arrange(sumz_wt_pval) %>%
    mutate(Rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(Rank, everything())

treatment_gene_set_pvals   <- read_feather(file.path(results_dir, "associations/treatment_response/gene_set/pvals.feather"))
treatment_gene_set_errors  <- read_feather(file.path(results_dir, "associations/treatment_response/gene_set/errors.feather"))
treatment_gene_set_effects <- read_feather(file.path(results_dir, "associations/treatment_response/gene_set/effects.feather"))

#--------------------------------
# x. Combined expr
#--------------------------------
gene_expr <- read_feather(file.path(results_dir, "expr/gene/expr-scaled-full.feather")) %>%
  column_to_rownames("symbol")

gene_set_expr <- read_feather(file.path(results_dir, "expr/gene_set/expr-scaled-full.feather")) %>%
  column_to_rownames("gene_set")

#preview_sample_ind <- sample(1:ncol(gene_expr), cfg$preview_num_samples)

# TEMP / use pipeline gene stats / "num_missing"
# num_nas <- apply(gene_expr, 1, function(x) { sum(is.na(x)) })
# mask <- num_nas <= cfg$preview_max_missing_samples
mask <- complete.cases(gene_expr)

preview_genes <- rownames(gene_expr)[mask]
preview_genes <- preview_genes[sample(1:length(preview_genes), cfg$preview_num_genes)]

#--------------------------------
# x. Co-expression
#--------------------------------
gene_coex <- read_feather(file.path(results_dir, "coex/gene/coex-scaled-5k.feather")) %>%
  column_to_rownames("symbol")

gene_set_coex <- read_feather(file.path(results_dir, "coex/gene_set/coex-scaled-5k.feather")) %>%
  column_to_rownames("gene_set")

gene_coex_opts <- gene_scores_all %>%
  filter(symbol %in% rownames(gene_coex)) %>%
  pull(symbol)

gene_set_coex_opts <- gene_set_scores_all %>%
  filter(gene_set %in% rownames(gene_set_coex)) %>%
  pull(gene_set)

#--------------------------------
# x. TCGA gene survival data
#--------------------------------
tcga <- read_feather("../scripts/data/hpa/hpa-cancer.feather")

#--------------------------------
# x. HMCL expr
#--------------------------------
expr_hmcl <- read_feather("../data/7.3/expr/expr_hmcl.feather")
hmcl_cells <- colnames(expr_hmcl)[-1]

log_info("Finished loading data..")
