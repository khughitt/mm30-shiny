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
library(plotly)
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(shinyWidgets)
library(survminer)
library(svglite)
library(tidyverse)
library(yaml)

source("R/plotting.R")

#############################
# Setup
#############################
options(spinner.color="#00bc8c")
options(scipen=2, digits=3)
set.seed(1)

enableBookmarking("url")

#log_threshold(DEBUG)
log_info("Initializing MM30 shiny app..")

tableOpts  <- list(pageLength=20, scrollX=TRUE)
selectOpts <- list(target="row", mode="single", selected=1)

cfg <- read_yaml("config/config-v7.3.yml")

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

pathway_mdata <- read_feather(file.path(results_dir, "metadata/pathways.feather"))

sample_mdata <- read_feather(file.path(results_dir, "metadata/samples.feather"))

covariate_mdata <- read_feather(file.path(results_dir, "metadata/covariates.feather"))

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

#--------------------------------
# 3. Overall survival
#--------------------------------
surv_os_gene_scores <- read_feather(file.path(results_dir, "scores/combined/survival_os/gene.feather")) %>%
    left_join(gene_mdata, by='symbol') %>%
    filter(num_missing <= cfg$max_missing$surv_os) %>%
    select(symbol, sumz_wt_pval, metafor_pval, description,
           chr_region, cell_cycle_phase, dgidb_categories) %>%
    arrange(sumz_wt_pval) %>%
    mutate(rank=dense_rank(sumz_wt_pval)) %>%
    select(rank, everything())

surv_os_gene_pvals   <- read_feather(file.path(results_dir, "associations/survival_os/gene/pvals.feather"))
surv_os_gene_errors  <- read_feather(file.path(results_dir, "associations/survival_os/gene/errors.feather"))
surv_os_gene_effects <- read_feather(file.path(results_dir, "associations/survival_os/gene/effects.feather"))

surv_os_pathway_scores <- read_feather(file.path(results_dir, "scores/combined/survival_os/pathway.feather")) %>%
    left_join(pathway_mdata, by='gene_set') %>%
    filter(num_missing <= cfg$max_missing$surv_os) %>%
    select(Pathway=gene_set, sumz_wt_pval, `P-value\n(metafor)`=metafor_pval, Collection=collection) %>%
    arrange(sumz_wt_pval) %>%
    mutate(rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(rank, everything())

surv_os_pathway_pvals   <- read_feather(file.path(results_dir, "associations/survival_os/pathway/pvals.feather"))
surv_os_pathway_errors  <- read_feather(file.path(results_dir, "associations/survival_os/pathway/errors.feather"))
surv_os_pathway_effects <- read_feather(file.path(results_dir, "associations/survival_os/pathway/effects.feather"))

#--------------------------------
# 3. Progression free survival 
#--------------------------------
gene_scores_surv_pfs <- read_feather(file.path(results_dir, "scores/combined/survival_pfs/gene.feather")) %>%
    left_join(gene_mdata, by='symbol') %>%
    filter(num_missing <= cfg$max_missing$surv_pfs) %>%
    select(Gene=symbol, sumz_wt_pval, `P-value\n(metafor)`=metafor_pval, Description=description,
           CHR=chr_region, `Cell Cycle`=cell_cycle_phase, `DGIdb\ncategories`=dgidb_categories) %>%
    arrange(sumz_wt_pval) %>%
    mutate(rank=dense_rank(sumz_wt_pval)) %>%
    rename(`P-value\n(metap)`=sumz_wt_pval) %>%
    select(rank, everything())


#--------------------------------
# x. Combined expr
#--------------------------------
expr_dat <- read_feather(file.path(results_dir, "expr/gene/expr_scaled.feather")) %>%
  column_to_rownames("symbol")

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
