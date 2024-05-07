#' Comparison of Powers for Different Tests under Different DBCD Randomization
#' Methods (Gaussian Responses)
#'
#' Compares the power of tests under different sample sizes for the same
#' treatment effects and design through matrices and plots.
#'
#' @usage power_comparison_Power_vs_n_GAUSSIAN(n_seq, nstart_seq, mu, sd,
#' nstop_mat, replication, group_allo, rho_func_index, rho_func, alpha,
#' sig_level)
#'
#' @param n_seq A sequence of numbers of patients.
#' The default is c(50, 100, 150, 200).
#' @param nstart_seq burn-in sample size of each arm. If NULL, n_seq/20 will be used.
#' @param mu A vector of mean response for each treatment arm
#' (where the first element refers to the control arm).
#' The length of mu should correspond to the number of arms.
#' The default is mu = c(4.5,5).
#' @param sd A vector of response standard deviations for each treatment arm.
#' (where the first element refers to the control arm).
#' The length of sd should correspond to the number of arms.
#' The default is sd = c(1.32, 0.72).
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
#' @param rho_func_index Supply a number of 1 or 2 indicting the
#' allocation function to use.
#' 1 = Zhang-Rosenberger allocation;
#' 2 (default) = Neyman allocation.
#' @param rho_func Supply a user-specified allocation function of Mean_RK and
#' SD_RK when rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. The default is 2.
#' @param sig_level Significant level. The default is 0.05.
#'
#'
#' @details 'power_comparison_Power_vs_n_GAUSSIAN' reads different sample sizes
#' as well as the corresponding burn-in size and outputs allocation, estimated
#' rates and powers.
#'
#' @return
#' \itemize{
#'   \item Allocation - Average and standard deviation (SD) of allocation distribution
#'   \item Estimation - Average and standard deviation of treatment effect
#'   \item Power_oneside - Average power of one-sided Welch T-test performed for each of the k-th arm against H0: mu_1>mu_k without multiplicity adjustment
#'   \item Power_aov - Average power of ANOVA test
#'   \item Plot - Four figures of results: 1) Allocation mean and SD, 2) Estimated mean response and SD, 3) Power of ANOVA test, 4) Power of one-sided Welch T-test
#' }
#'
#'
#' @export
#' @importFrom stats sd
#' @importFrom stats power
#' @importFrom stats lm
#' @import patchwork
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#'
#' ## Default setting
#' power_comparison_Power_vs_n_GAUSSIAN(
#' n_seq = seq(from = 50, to = 200, by = 50),
#' nstart_seq = round(seq(from = 50, to = 200, by = 50) / 20),
#' mu = c(4.5,5),
#' sd = c(1.32,0.72),
#' nstop_mat = NULL,
#' replication = 5,
#' group_allo = 1,
#' rho_func_index = 2,
#' rho_func = NULL,
#' alpha = 2,
#' sig_level = 0.05
#' )
#'
#'
power_comparison_Power_vs_n_GAUSSIAN = function(n_seq = seq(from = 50, to = 200, by = 50),
                                                nstart_seq = round(n_seq / 20),
                                                mu = c(4.5,5),
                                                sd = c(1.32,0.72),
                                                nstop_mat = NULL,
                                                replication = 100,
                                                group_allo = 1,
                                                rho_func_index = 2,
                                                rho_func = NULL,
                                                alpha = 2,
                                                sig_level = 0.05) {
  allocation <- estimate <- group <- samplesize <- mu_estimate <- power_aov <- NULL
  K = length(mu)
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
      simulation_main_GAUSSIAN(n = n,
                               nstart = nstart,
                               mu = mu,
                               sd = sd,
                               nstop = nstop,
                               replication = replication,
                               group_allo = group_allo,
                               rho_func_index = rho_func_index,
                               rho_func = rho_func,
                               alpha = alpha,
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
                      position = ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = estimate - sd, ymax = estimate + sd),
      width = dodge/5,
      position = ggplot2::position_dodge(dodge/1.1),
      color = "gray"
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

  ## plot mu-hat

  mu_est_data = outputmat[, colnames(outputmat) %in% c("samplesize", paste0("mu_estimate_mean", 1:K))]
  plotdata_mu_est = tidyr::gather(
    as.data.frame(mu_est_data),
    mu_estimate,
    estimate,
    paste0("mu_estimate_mean", 1:K),
    factor_key = TRUE
  )
  plotdata_mu_est = cbind(plotdata_mu_est, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_mu_est$mu_estimate
  ))) %>% dplyr::arrange(group,samplesize)

  mu_est_data_sd = outputmat[,colnames(outputmat) %in% c("samplesize", paste0("mu_estimate_sd",1:K))]
  plotdata_mu_est_sd = tidyr::gather(
    as.data.frame(mu_est_data_sd),
    mu_estimate,
    sd,
    paste0("mu_estimate_sd", 1:K),
    factor_key = TRUE
  )
  plotdata_mu_est_sd = cbind(plotdata_mu_est_sd, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_mu_est_sd$mu_estimate
  ))) %>% dplyr::arrange(group,samplesize)

  plotdata_mu_est_sd = plotdata_mu_est_sd[, c("samplesize", "sd", "group")]
  plotdata_mu_est = merge(plotdata_mu_est, plotdata_mu_est_sd, by = c("samplesize", "group")) %>% dplyr::arrange(group,samplesize)


  ggplot2::ggplot(data = plotdata_mu_est, ggplot2::aes(
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
      vjust = -0.55,
      color = "black",
      position = ggplot2::position_dodge(dodge/1.1),
      size = 3.5,
      angle = 30
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette = "Pastel1", name = "Trt") +
    ggplot2::ylim(min(0,min(plotdata_mu_est$estimate)-abs(min(plotdata_mu_est$estimate))/10),
                  max(0,max(plotdata_mu_est$estimate)+abs(max(plotdata_mu_est$estimate))/10)) +
    ggplot2::labs(title = "Estimated mean response",
                  x = "sample size",
                  y = latex2exp::TeX("$\\hat{\\mu}$")) +
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
      title = "Power - One sided Welch T-test",
      subtitle = "Reference to Trt 1 (placebo)",
      x = "sample size",
      y = "power"
    ) +
    ggplot2::scale_x_continuous(breaks = n_seq)  -> p_oneside

  ## plot ANOVA test
  plotdata_power_aov = as.data.frame(outputmat[, c("samplesize", "power_aov")])
  ggplot2::ggplot(data = plotdata_power_aov,
                  ggplot2::aes(x = samplesize, y = power_aov)) +
    ggplot2::geom_bar(
      stat = "identity",
      color = "black",
      position = ggplot2::position_dodge(),
      fill = "palegreen4"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(power_aov, 2)),
      vjust = -0.6,
      color = "black",
      position = ggplot2::position_dodge(dodge/1.3),
      size = 3.5,
      angle = 30
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0, 1.15) +
    ggplot2::labs(title = "Power - ANOVA test",
                  x = "sample size", y = "power") +
    ggplot2::scale_x_continuous(breaks = n_seq)  -> p_aov

  result = list(
    Allocation = plotdata_allo,
    # Allocation_plot = p_allocation,
    Estimation = plotdata_mu_est,
    Power_aov = plotdata_power_aov,
    # Power_aov_plot = p_aov
    Power_oneside = plotdata_power,
    # Power_oneside_plot = p_oneside,
    Plot = (p_allocation + p_est) / (p_aov + p_oneside)
  )
  return(result)
}


