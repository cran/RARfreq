#' Comparison of Powers for Sample Sizes under Different SEU Randomization
#' Methods (Binary Responses)
#'
#' Compares the power of tests under different sample sizes for the same
#' treatment effects and design through matrices and plots.
#'
#' @usage SEU_power_comparison_Power_vs_n(n_seq, nstart_seq, p, urn_comp,
#' nstop_mat, replication, group_allo, add_rule_index, add_rule, add_rule_full,
#' sig_level)
#'
#' @param n_seq A sequence of settings' number of patients.
#' The default is c(50, 100, 150, 200).
#' @param nstart_seq The burn-in sample size of each arm. The default is
#' n_seq/20 = c(2, 5, 8, 10).
#' @param p A vector of probabilities containing probabilities
#' for each treatment arm (where the first element refers to the control arm).
#' The length of p should correspond to the number of treatment arms.
#' The default is p = c(0.3, 0.3, 0.6).
#' @param urn_comp A vector of current urn composition. The default is NULL,
#' which indicates no ball in the urn.
#' @param nstop_mat A matrix of sample size stopping caps for each arm.
#' Each row corresponds to each n in n_seq, and each column represents each arm.
#' The trial stops if at least one arm reaches the corresponding cap.
#' The default is NULL, which means no cap.
#' @param replication the number of replications of the simulation. The default
#' is 100.
#' @param group_allo A number or a vector of group size(s) for allocation.
#' If a number is given, the allocation ratios will be updated for each batch of
#'  group_allo samples. If a vector is given, the allocation ratios will be
#'  updated sequentially in group according to the vector.
#'  The group_allo will be applied to all n (from each n_seq).
#'  Any value greater than n  will be omitted.
#' The default is group_allo=1, which is the same as group_allo = seq(nstart*length(p)+1,n).
#' @param add_rule_index Supply a number of 1, 2 or 3 indicting the
#' addition rules to target allocation functions.
#' 1 = randomized play-the-winner (RPW) rule that targets the urn allocation
#' 2 = the SEU model that targets Neyman allocation;
#' 3 = the SEU model that targets Rosenberger allocation;'
#' 4 = the SEU model that assigns probability of 0.6+1/K to winner at each step.
#' The default is 1.
#' @param add_rule Supply a user-specified addition rules function of x.df and
#' arms when add_rule_index is NULL. Default is NULL. (See SEU_BINARY_raw for
#' details on x.df and arms.)
#' @param add_rule_full Indicator of reference data for updating addition rule.
#' If TRUE, the addition rule is updated by full observation at each group
#' allocation. If FALSE,the addition rule is updated by each group observation.
#' The default is FALSE for add_rule_index=1 and TRUE otherwise.
#' @param sig_level Significant level (one-sided). The default is 0.05.
#'
#'
#' @details 'SEU_power_comparison_Power_vs_n' reads different sample sizes as
#' well as the corresponding burn-in size and outputs allocation, estimated
#' rates and powers.
#'
#' @return
#' \itemize{
#'   \item Allocation - Average and standard deviation (SD) of allocation distribution
#'   \item Estimation - Average and standard deviation of treatment effect
#'   \item Power_chisq - Average power of Chi-square test
#'   \item Power_oneside - Average power of one-sided Welch T-test performed for each of the k-th arm against H0: p_1>p_k without multiplicity adjustment
#'   \item Plot - Four figures of results: 1) Allocation mean and SD, 2) Estimated mean response and SD, 3) Power of Chi-square test, 4) Power of one-sided proportion test
#' }
#'
#' @export
#' @importFrom stats sd
#' @importFrom stats power
#' @import patchwork
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#'
#' ## Default setting
#' SEU_power_comparison_Power_vs_n(
#' n_seq = seq(from = 50, to = 200, by = 50),
#' nstart_seq = round(seq(from = 50, to = 200, by = 50) / 20),
#' p = c(0.3, 0.3, 0.6),
#' urn_comp = NULL,
#' nstop_mat = NULL,
#' replication = 4,
#' group_allo = 1,
#' add_rule_index = 1,
#' add_rule = NULL,
#' add_rule_full = FALSE,
#' sig_level = 0.05
#' )
#'
#'
SEU_power_comparison_Power_vs_n = function(n_seq = seq(from=50, to=200, by=50),
                                           nstart_seq = NULL,
                                           p = c(0.3, 0.3, 0.6),
                                           urn_comp = NULL,
                                           nstop_mat = NULL,
                                           replication = 100,
                                           group_allo = 1,
                                           add_rule_index = 1,
                                           add_rule = NULL,
                                           add_rule_full = NULL,
                                           sig_level = 0.05) {
  allocation <- estimate <- group <- samplesize <- p_estimate <- power_chisq <- NULL


  K = length(p)
  if(is.null(nstart_seq)) nstart_seq = round(n_seq/20) else
    if(length(nstart_seq)!=length(n_seq)) stop("The burn-in sample size is incampatible with n.")
  n_order = order(n_seq)
  n_seq = n_seq[n_order]
  nstart_seq = nstart_seq[n_order]
  if(!is.null(nstop_mat)) if(!all(dim(nstop_mat)==c(length(n_seq),K)))
    stop("The dimension of stopping cap matrix should be consistent with n_seq and p.")

  ## use simulation_main() to output power
  outputlist = list()
  outputmat = c()
  for (i_n in seq(n_seq)) {
    n = n_seq[i_n]
    nstart = nstart_seq[i_n]
    nstop = nstop_mat[i_n,]
    outputlist[[i_n]] = c(
      samplesize = n,
      SEU_simulation_main(n = n,
                          nstart = nstart,
                          p = p,
                          urn_comp = urn_comp,
                          nstop = nstop,
                          replication = replication,
                          group_allo = group_allo,
                          add_rule_index = add_rule_index,
                          add_rule = add_rule,
                          add_rule_full = add_rule_full,
                          sig_level = sig_level)
    )
    outputmat = rbind(outputmat, unlist(outputlist[[i_n]]))

  }

  ## plot allocation

  allocationmean = outputmat[, colnames(outputmat) %in% c("samplesize", paste0("allocation_mean", 1:K))]
  plotdata = tidyr::gather(
    as.data.frame(allocationmean),
    allocation,
    estimate,
    paste0("allocation_mean", 1:K),
    factor_key = TRUE
  )
  plotdata = cbind(plotdata, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata$allocation
  )))
  allocationsd = outputmat[, colnames(outputmat) %in% c("samplesize", paste0("allocation_sd", 1:K))]
  sddata = tidyr::gather(
    as.data.frame(allocationsd),
    allocation,
    sd,
    paste0("allocation_sd", 1:K),
    factor_key = TRUE
  )
  sddata = cbind(sddata, group = as.numeric(gsub(".+([0-9]+)", "\\1", sddata$allocation)))
  sddata = sddata[, c("samplesize", "sd", "group")]
  plotdata_allo = merge(plotdata, sddata, by = c("samplesize", "group")) %>% dplyr::arrange(group,samplesize)

  dodge = min(abs(diff(sort(n_seq))))

  ggplot2::ggplot(data = plotdata_allo, ggplot2::aes(
    x = samplesize,
    y = estimate,
    fill = factor(group)
  )) +
    ggplot2::geom_bar(stat = "identity",
                      color = "black",
                      position = ggplot2::position_dodge()
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = estimate - sd, ymax = estimate + sd),
      width = dodge/5,
      position = ggplot2::position_dodge(dodge/1.1),
      color = "dodgerblue3"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(estimate, 2)),
      vjust = -0.6,
      color = "black",
      position = ggplot2::position_dodge(dodge/1.1),
      size = 3.5,
      angle = 30
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette = "Paired", name = "Trt") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = "Allocation mean",
                  x = "sample size", y = "mean") +
    ggplot2::scale_x_continuous(breaks = n_seq)  -> p_allocation

  ## plot p-hat

  p_est_data = outputmat[, colnames(outputmat) %in% c("samplesize", paste0("p_estimate_mean", 1:K))]
  plotdata_p_est = tidyr::gather(
    as.data.frame(p_est_data),
    p_estimate,
    estimate,
    paste0("p_estimate_mean", 1:K),
    factor_key = TRUE
  )
  plotdata_p_est = cbind(plotdata_p_est, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_p_est$p_estimate
  ))) %>% dplyr::arrange(group,samplesize)


  p_est_data_sd = outputmat[, colnames(outputmat) %in% c("samplesize", paste0("p_estimate_sd", 1:K))]
  plotdata_p_est_sd = tidyr::gather(
    as.data.frame(p_est_data_sd),
    p_estimate,
    sd,
    paste0("p_estimate_sd", 1:K),
    factor_key = TRUE
  )
  plotdata_p_est_sd = cbind(plotdata_p_est_sd, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_p_est_sd$p_estimate
  ))) %>% dplyr::arrange(samplesize,group)

  plotdata_p_est_sd = plotdata_p_est_sd[, c("samplesize", "sd", "group")]
  plotdata_p_est = merge(plotdata_p_est, plotdata_p_est_sd, by = c("samplesize", "group")) %>% dplyr::arrange(group,samplesize)


  ggplot2::ggplot(data = plotdata_p_est, ggplot2::aes(
    x = samplesize,
    y = estimate,
    fill = factor(group)
  )) +
    ggplot2::geom_bar(stat = "identity",
                      color = "black",
                      position = ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-sd, ymax=estimate+sd),
                           width=dodge/5,
                           position=ggplot2::position_dodge(dodge/1.1),
                           color="gray") +
    ggplot2::geom_text(
      ggplot2::aes(label = round(estimate, 2)),
      vjust = -0.6,
      color = "black",
      position = ggplot2::position_dodge(dodge/1.1),
      size = 3.5,
      angle = 30
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette = "Pastel1", name = "Trt") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title = "Estimated success rate",
                  x = "sample size",
                  y = latex2exp::TeX("$\\hat{p}$")) +
    ggplot2::scale_x_continuous(breaks = n_seq)  -> p_est

  ## plot power of one-sided proportion test

  unadjustpower = outputmat[, colnames(outputmat) %in% c("samplesize", paste0("power_oneside.Trt ", 1:K))]
  plotdata = tidyr::gather(
    as.data.frame(unadjustpower),
    power,
    estimate,
    paste0("power_oneside.Trt ", 1:K),
    factor_key = TRUE
  )
  plotdata_power = cbind(plotdata, group = as.numeric(gsub(".+([0-9]+)", "\\1", plotdata$power))) %>%
    dplyr::arrange(group,samplesize)

  ggplot2::ggplot(data = plotdata_power, ggplot2::aes(
    x = samplesize,
    y = estimate,
    fill = factor(group)
  )) +
    ggplot2::geom_bar(
      data = plotdata_power,
      ggplot2::aes(
        x = samplesize,
        y = estimate,
        fill = factor(group)
      ),
      stat = "identity",
      color = "black",
      position = ggplot2::position_dodge()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(estimate, 2)),
      vjust = -0.6,
      color = "black",
      position = ggplot2::position_dodge(dodge/1.3),
      size = 3.5,
      angle = 30
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette = "GnBu", name = "Trt") +
    ggplot2::ylim(0, 1.15) +
    ggplot2::labs(
      title = "Power - One sided proportion test",
      subtitle = "Reference to Trt 1 (placebo)",
      x = "sample size",
      y = "power"
    ) +
    ggplot2::scale_x_continuous(breaks = n_seq)  -> p_oneside

  ## plot chi-sq test

  # ggplot() +
  #   geom_line(data=as.data.frame(outputmat[,c("samplesize","power_chisq")]),
  #             aes(x=samplesize, y=power_chisq), color='red')+
  #   theme_minimal() +
  #   ylim(0, 1) -> p_chisq
  plotdata_power_chisq = as.data.frame(outputmat[, c("samplesize", "power_chisq")])
  ggplot2::ggplot(data = plotdata_power_chisq,
                  ggplot2::aes(x = samplesize, y = power_chisq)) +
    ggplot2::geom_bar(
      stat = "identity",
      color = "black",
      position = ggplot2::position_dodge(),
      fill = "palegreen4"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(power_chisq, 2)),
      vjust = -0.6,
      color = "black",
      position = ggplot2::position_dodge(dodge/1.3),
      size = 3.5,
      angle = 30
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0, 1.15) +
    ggplot2::labs(title = "Power - Chi-square test",
                  x = "sample size", y = "power") +
    ggplot2::scale_x_continuous(breaks = n_seq)  -> p_chisq

  result = list(
    Allocation = plotdata_allo,
    # Allocation_plot = p_allocation,
    Power_chisq = plotdata_power_chisq,
    # Power_chisq_plot = p_chisq
    Power_oneside = plotdata_power,
    # Power_oneside_plot = p_oneside,
    Estimation = plotdata_p_est,
    Plot = (p_allocation + p_est) / (p_chisq + p_oneside)
  )
  return(result)
}


