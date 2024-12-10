################################################################################
#
# Plotting functions
#
################################################################################

#
# Creates a kaplan-meyer survival plot for specified upper/lower expression quantiles.
#
plot_survival <- function(dat, dataset, covariate, feat_name, time_units, expr_cutoffs, color_pal) {
  # divide expression into quantiles
  cutoff <- as.numeric(expr_cutoffs)

  expr_quantiles <- quantile(dat$feature, c(cutoff / 100, 1 - (cutoff / 100)), na.rm=TRUE)

  dat$Expression <- ""
  dat$Expression[dat$feature <= expr_quantiles[1]] <- paste0("Lower ", names(expr_quantiles)[1])
  dat$Expression[dat$feature >= expr_quantiles[2]] <- paste0("Upper ", names(expr_quantiles)[1])

  # determine units to use
  # if (dataset %in% c("MMRF", "GSE7039", "GSE57317", "GSE9782")) {
  #   time_units <- "days"
  # } else if (dataset %in% c("GSE24080")) {
  #   time_units <- "weeks"
  # } else if (dataset %in% c("GSE19784")) {
  #   time_units <- "months"
  # } else {
  #   time_units <- "?"
  # }

  # drop all data except for upper and lower quantiles
  dat <- dat %>%
    filter(Expression != "")

  num_samples <- nrow(dat)

  dat$Expression <- factor(dat$Expression)

  cov_label <- str_to_title(gsub("_", " ", covariate))
  plt_title <- sprintf("%s: %s vs. %s Expression (n=%d)", dataset, cov_label, feat_name, num_samples)

  # perform fit on binarized data
  fit <- survival::survfit(survival::Surv(time, event) ~ Expression, data=dat)

  # display a kaplan meier plot for result
  survminer::ggsurvplot(fit, data=dat, ggtheme=theme_pubr(base_size=16), palette=color_pal,
                        title=plt_title, xlab=sprintf("Time (%s)", time_units),
                        legend="bottom", legend.title="Legend")

}

# Creates a violin plot for a given categorical response variable.
plot_categorical <- function(dat, dataset, covariate, feat_name, color_pal) {
  # drop any entries with missing values
  dat <- dat[!is.na(dat$response), ]

  # draw violin + jitter plot
  set.seed(1)

  ncat <- nlevels(dat$response)

  # if more factor levels exist than colors, expand palette
  if (ncat > length(color_pal)) {
    color_pal <- colorRampPalette(color_pal)(nlevels(dat$response))
  }

  plt <- ggplot(dat, aes(x=response, y=feature)) +
    geom_violin(aes(fill=response, color=response), alpha=0.5, draw_quantiles=c(0.5)) +
    geom_boxplot(width=0.1) +
    geom_jitter(aes(color=response), alpha=0.8) +
    scale_fill_manual(values=color_pal) +
    scale_color_manual(values=color_pal, guide="none") +
    theme_pubr(base_size=16) +
    xlab(covariate) +
    ylab(sprintf("%s expression", feat_name)) +
    ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, covariate))

  # work-around to hide part of legend with plotly/ggplotly
  # https://github.com/plotly/plotly.R/issues/572#issuecomment-876540008
  plt %>%
    style(plt, showlegend=FALSE, traces=(ncat + 1):(ncat * 2 + 1))
}

# Creates a scatter plot for differential expression covariates
plot_deseq <- function(dat, dataset, covariate, feat_name, color_pal) {
  # drop any entries with missing values
  #dat <- dat[!is.na(dat$response), ]

  # if more factor levels exist than colors, expand palette
  if (nlevels(dat$response) > length(color_pal)) {
    color_pal <- colorRampPalette(color_pal)(nlevels(dat$response))
  }

  var1 <- colnames(dat)[2]

  # single variable models
  if (ncol(dat) == 2) {
    ggplot(dat) +
      geom_col(aes(x=.data[[var1]], y=feature, fill=.data[[var1]], color=.data[[var1]])) +
      scale_fill_manual(values=color_pal) +
      scale_color_manual(values=color_pal) +
      ggtitle(sprintf("%s: %s vs. %s", dataset, feat_name, var1)) +
      theme_pubr(base_size=16) +
      xlab(covariate) +
      ylab(sprintf("%s expression", feat_name))
  }
}
