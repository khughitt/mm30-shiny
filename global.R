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

data_dir <- "../data"

scores_dir <- file.path(data_dir, "scores")
expr_dir <- file.path(data_dir, "expr")
coex_dir <- file.path(data_dir, "coex")
hmcl_dir <- file.path(data_dir, "hmcl")
hpa_dir <- file.path(data_dir, "hpa")
assoc_dir <- file.path(data_dir, "assoc")
other_dir <- file.path(data_dir, "other")
disease_stage_dir <- file.path(data_dir, "disease_stage")
transitions_dir <- file.path(data_dir, "transitions")
metadata_dir <- file.path(data_dir, "metadata")

#--------------------------
# 1. Metadata
#--------------------------
gene_mdata <- read_feather(file.path(metadata_dir, "genes.feather"))
gene_set_mdata <- read_feather(file.path(metadata_dir, "gene_sets.feather"))
sample_mdata <- read_feather(file.path(metadata_dir, "samples.feather"))
covariate_mdata <- read_feather(file.path(metadata_dir, "covariates.feather"))
dataset_mdata <- read_feather(file.path(metadata_dir, "datasets.feather"))
indiv_mdata <- readRDS(file.path(metadata_dir, "indiv_datasets.rds"))

surv_os_datasets <- covariate_mdata %>%
  filter(phenotype == "overall_survival") %>%
  pull(dataset) %>%
  sort(TRUE)

surv_pfs_datasets <- covariate_mdata %>%
  filter(phenotype == "prog_free_survival") %>%
  pull(dataset) %>%
  sort(TRUE)

stage_datasets <- covariate_mdata %>%
  filter(phenotype=="disease_stage") %>%
  pull(dataset) %>%
  sort(TRUE)

treatment_datasets <- covariate_mdata %>%
  filter(phenotype=="treatment_response") %>%
  pull(dataset) %>%
  sort(TRUE)

#--------------------------------
# 2. "All"
#--------------------------------
gene_scores_all <- read_feather(file.path(scores_dir, "all/gene.feather"))
gene_set_scores_all <- read_feather(file.path(scores_dir, "all/gene_set.feather"))

#--------------------------------
# 3. Overall survival
#--------------------------------
surv_os_gene_scores <- read_feather(file.path(scores_dir, "survival_os/gene.feather"))
surv_os_gene_pvals   <- read_feather(file.path(assoc_dir, "survival_os/gene/pvals.feather"))
surv_os_gene_errors  <- read_feather(file.path(assoc_dir, "survival_os/gene/errors.feather"))
surv_os_gene_effects <- read_feather(file.path(assoc_dir, "survival_os/gene/effects.feather"))

surv_os_gene_set_scores <- read_feather(file.path(scores_dir, "survival_os/gene_set.feather"))
surv_os_gene_set_pvals   <- read_feather(file.path(assoc_dir, "survival_os/gene_set/pvals.feather"))
surv_os_gene_set_errors  <- read_feather(file.path(assoc_dir, "survival_os/gene_set/errors.feather"))
surv_os_gene_set_effects <- read_feather(file.path(assoc_dir, "survival_os/gene_set/effects.feather"))

#--------------------------------
# 4. Progression free survival
#--------------------------------
surv_pfs_gene_scores <- read_feather(file.path(scores_dir, "survival_pfs/gene.feather"))
surv_pfs_gene_pvals   <- read_feather(file.path(assoc_dir, "survival_pfs/gene/pvals.feather"))
surv_pfs_gene_errors  <- read_feather(file.path(assoc_dir, "survival_pfs/gene/errors.feather"))
surv_pfs_gene_effects <- read_feather(file.path(assoc_dir, "survival_pfs/gene/effects.feather"))

surv_pfs_gene_set_scores <- read_feather(file.path(scores_dir, "survival_pfs/gene_set.feather"))
surv_pfs_gene_set_pvals   <- read_feather(file.path(assoc_dir, "survival_pfs/gene_set/pvals.feather"))
surv_pfs_gene_set_errors  <- read_feather(file.path(assoc_dir, "survival_pfs/gene_set/errors.feather"))
surv_pfs_gene_set_effects <- read_feather(file.path(assoc_dir, "survival_pfs/gene_set/effects.feather"))

#--------------------------------
# x. Disease stage
#--------------------------------
gene_stage_scores <- read_feather(file.path(scores_dir, "disease_stage/gene.feather"))
gene_stage_pvals   <- read_feather(file.path(assoc_dir, "disease_stage/gene/pvals.feather"))
gene_stage_errors  <- read_feather(file.path(assoc_dir, "disease_stage/gene/errors.feather"))
gene_stage_effects <- read_feather(file.path(assoc_dir, "disease_stage/gene/effects.feather"))

gene_stage_scaled_expr <- read_feather(file.path(expr_dir, "gene", "disease_stage.feather"))

gene_set_stage_scores <- read_feather(file.path(scores_dir, "disease_stage/gene_set.feather"))
gene_set_stage_pvals   <- read_feather(file.path(assoc_dir, "disease_stage/gene_set/pvals.feather"))
gene_set_stage_errors  <- read_feather(file.path(assoc_dir, "disease_stage/gene_set/errors.feather"))
gene_set_stage_effects <- read_feather(file.path(assoc_dir, "disease_stage/gene_set/effects.feather"))

gene_set_stage_scaled_expr <- read_feather(file.path(expr_dir, "gene_set", "disease_stage.feather"))

gene_transitions <- readRDS(file.path(transitions_dir, "gene.rds"))
gene_set_transitions <- readRDS(file.path(transitions_dir, "gene_set.rds"))
transitions <- names(gene_transitions)

# number of datasets for which each of the transitions was observed
trans_counts <- lapply(gene_transitions, function(x) {
  # exclude "symbol", "num_datasets", "mean_change" and "median_change" columns from counts
  ncol(x) - 4
})

transition_select_opts <- c(
  "Healthy_MGUS",
  "MGUS_SMM",
  "SMM_MM",
  "MM_RRMM",
  "early_vs_late",
  "before_vs_after_smm",
  "before_vs_after_relapse"
)
transition_select_labels <- c(sprintf("Healthy vs. MGUS (n=%d)", trans_counts[["Healthy_MGUS"]]),
                              sprintf("MGUS vs. SMM (n=%d)", trans_counts[["MGUS_SMM"]]),
                              sprintf("SMM vs. MM (n=%d)", trans_counts[["SMM_MM"]]),
                              sprintf("MM vs. RRMM (n=%d)", trans_counts[["MM_RRMM"]]),
                              sprintf("Healthy/MGUS/SMM vs. MM/RRMM (n=%d)", trans_counts[["early_vs_late"]]),
                              sprintf("Healthy/MGUS vs. SMM/MM/RRMM (n=%d)", trans_counts[["before_vs_after_smm"]]),
                              sprintf("Pre- vs. Post-relapse (n=%d)", trans_counts[["before_vs_after_relapse"]]))

names(transition_select_opts) <- transition_select_labels

# mapping from transition names to stages
transition_stages <- list(
  "Healthy_MGUS"=c("Healthy", "MGUS"),
  "MGUS_SMM"=c("MGUS", "SMM"),
  "SMM_MM"=c("SMM", "MM"),
  "MM_RRMM"=c("MM", "RRMM"),
  "early_vs_late"=c("Healthy", "MGUS", "SMM", "MM", "RRMM"),
  "before_vs_after_smm"=c("Healthy", "MGUS", "SMM", "MM", "RRMM"),
  "before_vs_after_relapse"=c("Healthy", "MGUS", "SMM", "MM", "RRMM")
)

#--------------------------------
# x. Treatment
#--------------------------------
treatment_gene_scores <- read_feather(file.path(scores_dir, "treatment_response/gene.feather"))
treatment_gene_pvals   <- read_feather(file.path(assoc_dir, "treatment_response/gene/pvals.feather"))
treatment_gene_errors  <- read_feather(file.path(assoc_dir, "treatment_response/gene/errors.feather"))
treatment_gene_effects <- read_feather(file.path(assoc_dir, "treatment_response/gene/effects.feather"))

treatment_gene_set_scores <- read_feather(file.path(scores_dir, "treatment_response/gene_set.feather"))
treatment_gene_set_pvals   <- read_feather(file.path(assoc_dir, "treatment_response/gene_set/pvals.feather"))
treatment_gene_set_errors  <- read_feather(file.path(assoc_dir, "treatment_response/gene_set/errors.feather"))
treatment_gene_set_effects <- read_feather(file.path(assoc_dir, "treatment_response/gene_set/effects.feather"))

#--------------------------------
# x. Combined expr
#--------------------------------
gene_expr <- read_feather(file.path(expr_dir, "gene/scaled.feather")) %>%
  column_to_rownames("symbol")

gene_set_expr <- read_feather(file.path(expr_dir, "gene_set/scaled.feather")) %>%
  column_to_rownames("gene_set")

preview_genes <- readRDS(file.path(other_dir, "preview_genes.rds"))

#--------------------------------
# x. Co-expression (?)
#--------------------------------
# gene_coex <- read_feather(file.path(coex_dir, "gene/coex-scaled-5k.feather")) %>%
#   column_to_rownames("symbol")

# gene_coex_opts <- gene_scores_all %>%
#   filter(symbol %in% rownames(gene_coex)) %>%
#   pull(symbol)
gene_coex_opts <- gene_scores_all %>%
  pull(symbol)

# gene_set_coex_opts <- gene_set_scores_all %>%
#   filter(gene_set %in% rownames(gene_set_coex)) %>%
#   pull(gene_set)
gene_set_coex_opts <- gene_set_scores_all %>%
  pull(gene_set)

#--------------------------------
# x. TCGA gene survival data
#--------------------------------
tcga <- read_feather(file.path(hpa_dir, "hpa-cancer.feather"))

#--------------------------------
# x. HMCL expr
#--------------------------------
expr_hmcl <- read_feather(file.path(hmcl_dir, "expr.feather"))
hmcl_cells <- colnames(expr_hmcl)[-1]

log_info("Finished loading data..")
