#' Comparison of Powers for Treatment Effects under Different DBCD Randomization
#' methods (Gaussian Responses)
#'
#' Compares the power of a 2-arm design under different treatment effects
#' for the same sample size and placebo effect through matrices and plots.
#'
#' @usage power_comparison_Power_vs_Trt_GAUSSIAN(n, nstart, mu_pbo, sd_pbo,
#' mu_trt, sd_trt, nstop, replication, group_allo, rho_func_index, rho_func,
#' alpha, sig_level)
#'
#' @param n The number of patients. Default is 100.
#' @param nstart Burn-in sample size of each arm. Default is n/20.
#' @param mu_pbo Mean response of placebo (control) arm. Default is 4.5.
#' @param sd_pbo Response standard deviation of placebo (control) arm. Default
#' is 1.32.
#' @param mu_trt A vector containing mean responses of treatment.
#' @param sd_trt A vector containing response standard deviations of treatment.
#' @param nstop A vector of stopping cap of sample size for each arm. The trial
#' stops if at least one arm reaches the corresponding cap.
#' The default is NULL, which means no cap.
#' @param replication the number of replications of the simulation. The default
#' is 100.
#' @param group_allo A number or a vector of group size(s) for allocation.
#' If a number is given, the allocation ratios will be updated for each batch of
#'  group_allo samples. If a vector is given, the allocation ratios will be
#'  updated sequentially in group according to the vector.
#'  Any value greater than n will be omitted.
#' The default is group_allo=1, which is the same as group_allo = seq(nstart*length(p)+1,n).
#' @param rho_func_index Supply a number of 1 or 2 indicting the
#' allocation function to use.
#' 1 = Zhang-Rosenberger allocation (2-arm allocation only);
#' 2 (default) = Neyman allocation.
#' @param rho_func Supply a user-specified allocation function of Mean_RK and
#' SD_RK when rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. Default is 2.
#' @param sig_level Significant level. Default is 0.05.
#'
#'
#' @details 'power_comparison_Power_vs_Trt_GAUSSIAN' reads different treatment
#' effects and outputs allocation, estimated rates and powers.
#'
#' @return
#' \itemize{
#'   \item Allocation - Average and standard deviation (SD) of allocation distribution
#'   \item Estimation - Average and standard deviation of treatment effect
#'   \item Power - Average power: 1) ANOVA test, 2) one-sided Welch T-test performed for each of the k-th arm against H0: mu_1>mu_k without multiplicity adjustment
#'   \item Plot - Three figures of results: 1) Allocation mean and SD, 2) Estimated mean response and SD, 3) Power of ANOVA test and power of one-sided Welch T-test
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
#' ## Default setting
#' power_comparison_Power_vs_Trt_GAUSSIAN(
#' n = 100,
#' nstart = round(100/20),
#' mu_pbo = 4.5,
#' sd_pbo = 1.32,
#' mu_trt = seq(0,0.6,0.2)+4.5,
#' sd_trt = rep(1.32,4),
#' nstop = NULL,
#' replication = 5,
#' group_allo = 1,
#' rho_func_index = 2,
#' rho_func = NULL,
#' alpha = 2,
#' sig_level = 0.05
#' )
#'
#'
power_comparison_Power_vs_Trt_GAUSSIAN = function(n = 100,
                                                  nstart = NULL,
                                                  mu_pbo = 4.5,
                                                  sd_pbo = 1.32,
                                                  mu_trt = seq(0,0.6,0.2)+4.5,
                                                  sd_trt = rep(1.32,4),
                                                  nstop = NULL,
                                                  replication = 100,
                                                  group_allo = 1,
                                                  rho_func_index = 2,
                                                  rho_func = NULL,
                                                  alpha = 2,
                                                  sig_level = 0.05){
  allocation <- estimate <- group <- mu_estimate <- test <- NULL

  K = 2
  # delta=sort(delta)
  # mu_trt = mu_pbo + delta
  # sd_trt = sd_pbo * sd_factor
  mu_trt = sort(mu_trt)
  delta = mu_trt - mu_pbo
  if(any(sd_trt<0) ) stop(
    "Treatment response sdandard deviation 'sd_trt' should be positive!")
  if(length(unique(delta)) != length(delta)) stop("Please assign unique 'mu_trt'.")
  if(length(delta)!=length(sd_trt)) stop(
    "The mean responses of treatment 'mu_trt' should share the same lengths as that of the response standard deviations 'sd_trt'."
  )
  if(is.null(nstart)) nstart = round(n/20)
  nstart = round(nstart)


  ## use simulation_main() to output power
  outputlist = list()
  outputmat = c()
  for(i_trt in seq(mu_trt)){
    mu = c(mu_pbo,mu_trt[i_trt])
    sd = c(sd_pbo,sd_trt[i_trt])
    outputlist[[i_trt]] = c(mu_trt = mu_trt[i_trt], sd_trt = sd_trt[i_trt],
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
                                                     sig_level = sig_level))
    outputmat = rbind(outputmat,unlist(outputlist[[i_trt]]))
  }

  ## plot allocation

  allocationmean = outputmat[,colnames(outputmat) %in% c("mu_trt", "sd_trt", paste0("allocation_mean",1:K))]
  plotdata = tidyr::gather(as.data.frame(allocationmean), allocation, estimate, paste0("allocation_mean",1:K),  factor_key=TRUE)
  plotdata = cbind(plotdata, group = as.numeric(gsub(".+([0-9]+)", "\\1", plotdata$allocation)))
  allocationsd = outputmat[,colnames(outputmat) %in% c("mu_trt", "sd_trt", paste0("allocation_sd",1:K))]
  sddata = tidyr::gather(as.data.frame(allocationsd), allocation, sd, paste0("allocation_sd",1:K),  factor_key=TRUE)
  sddata = cbind(sddata, group = as.numeric(gsub(".+([0-9]+)", "\\1", sddata$allocation)))
  sddata = sddata[,c("mu_trt", "sd_trt","sd","group")]
  plotdata_allo = merge(plotdata,sddata,by=c("mu_trt", "sd_trt","group")) %>% arrange(group)

  dodge = min(abs(diff(sort(delta))))

  ggplot2::ggplot(data=plotdata_allo, ggplot2::aes(x=mu_trt, y=estimate, fill=factor(group))) +
    ggplot2::geom_bar(stat="identity",color="black", position=ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-sd, ymax=estimate+sd), width=dodge/5,
                           position=ggplot2::position_dodge(dodge/1.1), color="gray") +
    ggplot2::geom_text(ggplot2::aes(label=round(estimate,2)), vjust=-0.6, color="black",
                       position = ggplot2::position_dodge(dodge/1.1), size=3.5, angle=30) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Paired", name="Trt") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title="Allocation mean",
                  subtitle = paste0("mu_pbo = ",mu_pbo, "\nsd_trt = [", paste(sd_trt,collapse = ","), "]"),
                  x ="mean response of trt", y = "mean") +
    ggplot2::scale_x_continuous(breaks = mu_trt)  -> p_allocation

  ## plot mu-hat

  mu_est_data = outputmat[,colnames(outputmat) %in% c("mu_trt", "sd_trt", paste0("mu_estimate_mean",1:K))]
  plotdata_mu_est = tidyr::gather(as.data.frame(mu_est_data), mu_estimate, estimate, paste0("mu_estimate_mean",1:K),  factor_key=TRUE)
  plotdata_mu_est = cbind(plotdata_mu_est, group = as.numeric(gsub(".+([0-9]+)", "\\1", plotdata_mu_est$mu_estimate)))

  mu_est_data_sd = outputmat[,colnames(outputmat) %in% c("mu_trt", "sd_trt", paste0("mu_estimate_sd",1:K))]
  plotdata_mu_est_sd = tidyr::gather(
    as.data.frame(mu_est_data_sd),
    mu_estimate,
    sd,
    paste0("mu_estimate_sd", 1:K),
    factor_key = TRUE
  )
  plotdata_mu_est_sd = cbind(plotdata_mu_est_sd, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_mu_est_sd$mu_estimate
  ))) %>% dplyr::arrange(group)

  plotdata_mu_est_sd = plotdata_mu_est_sd[, c("mu_trt", "sd_trt", "sd", "group")]
  plotdata_mu_est = merge(plotdata_mu_est, plotdata_mu_est_sd, by = c("mu_trt", "sd_trt", "group")) %>% dplyr::arrange(group)

  ggplot2::ggplot(data=plotdata_mu_est, ggplot2::aes(x=mu_trt, y=estimate, fill=factor(group))) +
    ggplot2::geom_bar(stat="identity",color="black", position=ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-sd, ymax=estimate+sd), width=dodge/5,
                           position=ggplot2::position_dodge(dodge/1.1), color="gray") +
    ggplot2::geom_text(ggplot2::aes(label=round(estimate,2)), vjust=-0.6, color="black",
                       position = ggplot2::position_dodge(dodge/1.1), size=3.5, angle=30) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Pastel1", name="Trt") +
    ggplot2::ylim(min(0,min(plotdata_mu_est$estimate)-abs(min(plotdata_mu_est$estimate))/10),
                  max(0,max(plotdata_mu_est$estimate)+abs(max(plotdata_mu_est$estimate))/10)) +
    ggplot2::labs(title="Estimated mean response",
                  subtitle = paste0("mu_pbo = ",mu_pbo, "\nsd_trt = [", paste(sd_trt,collapse = ","), "]"),
                  x ="mean response of trt", y = latex2exp::TeX("$\\hat{\\mu}$")) +
    ggplot2::scale_x_continuous(breaks = mu_trt) -> p_est

  ## plot two tests in one plot, x-axis = difference of mean trt effect
  powers = as.data.frame(outputmat[,c("mu_trt", "sd_trt","power_aov","power_oneside.Trt 2")])
  colnames(powers)[3] <- "ANOVA test"
  colnames(powers)[4] <- "One-sided Welch T-test"
  plotdata_power = tidyr::gather(powers, test, power,c("ANOVA test","One-sided Welch T-test"),  factor_key=TRUE)

  ggplot2::ggplot(data= plotdata_power,
                  ggplot2::aes(x=mu_trt, y=power, color=test)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() + ggplot2::geom_point() +
    ggplot2::scale_color_brewer(palette="Dark2", name="Tests") +
    ggplot2::ylim(0, 1.) +
    ggplot2::labs(title="Power v.s Treatment response mean",
                  subtitle= paste0("Reference to Trt 1 (placebo) ", "mu_pbo = ",mu_pbo, "\nsd_trt = [", paste(sd_trt,collapse = ","), "]"),
                  x ="mean response of trt", y = "power") +
    ggplot2::theme(legend.position="top",
                   legend.justification="right",
                   legend.margin=ggplot2::margin(0,0,0,0),
                   legend.box.margin=ggplot2::margin(-10,10,-10,-10)) +
    ggplot2::scale_x_continuous(breaks = mu_trt)  -> combined_power_plot

  result = list(
    Allocation = plotdata_allo,
    Estimation = plotdata_mu_est,
    Power = plotdata_power,
    Plot = (p_allocation + p_est) / combined_power_plot
  )
  return(result)
}

