#
# creates a list of sample metadata dataframes for each dataset
#
get_sample_metadata_list <- function(mmrf_dir, geo_dir) {
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

  return(indiv_mdata)
}
