#' Comparison of Powers for Treatment Effects under Different DBCD Randomization
#' Methods (Binary Responses)
#'
#' Compares the power of a 2-arm design under different treatment effects
#' for the same sample size and placebo effect through matrices and plots.
#'
#' @usage power_comparison_Power_vs_Trt(n, nstart, p_pbo, p_trt, nstop,
#' replication, group_allo, rho_func_index, rho_func, alpha, sig_level)
#'
#' @param n The number of patients. The default is 100.
#' @param nstart Burn-in sample size of each arm. The default is n/20.
#' @param p_pbo Success rate of placebo (control) arm. The default is 0.3.
#' @param p_trt A vector containing success rates of treatment arm.
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
#' @param rho_func_index Supply a number of 1, 2 or 3 indicting the
#' allocation function to use.
#' 1 = Wei's allocation
#' 2 = Neyman allocation;
#' 3 (default) = Rosenberger allocation;
#' @param rho_func Supply a user-specified allocation function of Mean_RK and
#' SD_RK when rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. The default is 2.
#' @param sig_level Significant level. The default is 0.05.
#'
#'
#' @details 'power_comparison_Power_vs_Trt' reads different treatment effects
#' and outputs allocation, estimated rates and powers.
#'
#' @return
#' \itemize{
#'   \item Allocation - Average and standard deviation (SD) of allocation distribution
#'   \item Estimation - Average and standard deviation of treatment effect
#'   \item Power - Average power: 1) Chi-square test, 2) one-sided proportion test performed for each of the k-th arm against H0: p_1>p_k without multiplicity adjustment
#'   \item Plot - Three figures of results: 1) Allocation mean and SD, 2) Estimated mean response and SD, 3) Power of Chi-square test and power of one-sided proportion test
#' }
#'
#'
#' @export
#' @importFrom stats sd
#' @importFrom stats power
#' @import patchwork
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#' ## Default setting
#' power_comparison_Power_vs_Trt(
#' n = 100,
#' nstart = round(100/20),
#' p_pbo = 0.3,
#' p_trt = seq(0,0.3,0.1)+0.3,
#' nstop = NULL,
#' replication = 5,
#' group_allo = 1,
#' rho_func_index = 3,
#' rho_func = NULL,
#' alpha = 2,
#' sig_level = 0.05
#' )
#'
#'
power_comparison_Power_vs_Trt = function(n = 100,
                                         nstart = NULL,
                                         p_pbo = 0.3,
                                         p_trt = seq(0,0.3,0.1)+0.3,
                                         nstop = NULL,
                                         replication = 100,
                                         group_allo = 1,
                                         rho_func_index = 3,
                                         rho_func = NULL,
                                         alpha = 2,
                                         sig_level = 0.05){
  allocation <- estimate <- group <- p_estimate <- test <- NULL

  K = 2
  # delta = sort(delta)
  # p_trt = p_pbo + delta
  p_trt = sort(p_trt)
  delta = p_trt - p_pbo
  if(any(p_trt>1) | any(p_trt<0) | p_pbo<0 | p_pbo>1) stop(
    "Please check placebo response rate 'p_pbo' and reatment response rate 'p_trt'!")
  if(is.null(nstart)) nstart = round(n/20)
  nstart = round(nstart)

  ## use simulation_main() to output power
  outputlist = list()
  outputmat = c()
  for(i_p in seq(p_trt)){
    p = c(p_pbo,p_trt[i_p])
    outputlist[[i_p]] = c(p_trt = p_trt[i_p],
                          simulation_main(n = n,
                                          nstart = nstart,
                                          p = p,
                                          nstop = nstop,
                                          replication = replication,
                                          group_allo = group_allo,
                                          rho_func_index = rho_func_index,
                                          rho_func = rho_func,
                                          alpha = alpha,
                                          sig_level = sig_level))
    outputmat = rbind(outputmat,unlist(outputlist[[i_p]]))
  }

  ## plot allocation

  allocationmean = outputmat[,colnames(outputmat) %in% c("p_trt", paste0("allocation_mean",1:K))]
  plotdata = tidyr::gather(as.data.frame(allocationmean), allocation, estimate, paste0("allocation_mean",1:K),  factor_key=TRUE)
  plotdata = cbind(plotdata, group = as.numeric(gsub(".+([0-9]+)", "\\1", plotdata$allocation)))
  allocationsd = outputmat[,colnames(outputmat) %in% c("p_trt", paste0("allocation_sd",1:K))]
  sddata = tidyr::gather(as.data.frame(allocationsd), allocation, sd, paste0("allocation_sd",1:K),  factor_key=TRUE)
  sddata = cbind(sddata, group = as.numeric(gsub(".+([0-9]+)", "\\1", sddata$allocation)))
  sddata = sddata[,c("p_trt","sd","group")]
  plotdata_allo = merge(plotdata,sddata,by=c("p_trt","group")) %>% arrange(group)

  dodge = min(abs(diff(sort(delta))))

  ggplot2::ggplot(data=plotdata_allo, ggplot2::aes(x=p_trt, y=estimate, fill=factor(group))) +
    ggplot2::geom_bar(stat="identity",color="black", position=ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-sd, ymax=estimate+sd), width=dodge/5,
                           position=ggplot2::position_dodge(dodge/1.1), color="dodgerblue3") +
    ggplot2::geom_text(ggplot2::aes(label=round(estimate,2)), vjust=-0.6, color="black",
              position = ggplot2::position_dodge(dodge/1.1), size=3.5, angle=30) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Paired", name="Trt") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title="Allocation mean",
         subtitle = paste0("p_pbo = ",p_pbo),
         x ="success rate of trt", y = "mean") +
    ggplot2::scale_x_continuous(breaks = p_trt)  -> p_allocation

  ## plot p-hat
#
#   p_est_data = outputmat[,colnames(outputmat) %in% c("p_trt", paste0("p_estimate_mean",1:K))]
#   plotdata_p_est = tidyr::gather(as.data.frame(p_est_data), p_estimate, estimate, paste0("p_estimate_mean",1:K),  factor_key=TRUE)
#   plotdata_p_est = cbind(plotdata_p_est, group = as.numeric(gsub(".+([0-9]+)", "\\1", plotdata_p_est$p_estimate)))

  p_est_data = outputmat[, colnames(outputmat) %in% c("p_trt", paste0("p_estimate_mean", 1:K))]
  plotdata_p_est = tidyr::gather(
    as.data.frame(p_est_data),
    p_estimate,
    estimate,
    paste0("p_estimate_mean", 1:K),
    factor_key = TRUE
  )
  plotdata_p_est = cbind(plotdata_p_est, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_p_est$p_estimate
  ))) %>% dplyr::arrange(group)


  p_est_data_sd = outputmat[, colnames(outputmat) %in% c("p_trt", paste0("p_estimate_sd", 1:K))]
  plotdata_p_est_sd = tidyr::gather(
    as.data.frame(p_est_data_sd),
    p_estimate,
    sd,
    paste0("p_estimate_sd", 1:K),
    factor_key = TRUE
  )
  plotdata_p_est_sd = cbind(plotdata_p_est_sd, group = as.numeric(gsub(
    ".+([0-9]+)", "\\1", plotdata_p_est_sd$p_estimate
  ))) %>% dplyr::arrange(group)

  plotdata_p_est_sd = plotdata_p_est_sd[, c("p_trt", "sd", "group")]
  plotdata_p_est = merge(plotdata_p_est, plotdata_p_est_sd, by = c("p_trt", "group")) %>% dplyr::arrange(group)
  ###

  ggplot2::ggplot(data=plotdata_p_est, ggplot2::aes(x=p_trt, y=estimate, fill=factor(group))) +
    ggplot2::geom_bar(stat="identity",color="black", position=ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-sd, ymax=estimate+sd), width=dodge/5,
                           position=ggplot2::position_dodge(dodge/1.1), color="gray") +
    ggplot2::geom_text(ggplot2::aes(label=round(estimate,2)), vjust=-0.6, color="black",
              position = ggplot2::position_dodge(dodge/1.1), size=3.5, angle=30) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer(palette="Pastel1", name="Trt") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(title="Estimated success rate",
         subtitle = paste0("p_pbo = ",p_pbo),
         x ="success rate of trt", y = latex2exp::TeX("$\\hat{p}$")) +
    ggplot2::scale_x_continuous(breaks = p_trt) -> p_est

  # ## plot power of one-sided proportion test
  #
  # unadjustpower = outputmat[,colnames(outputmat) %in% c("p_trt", paste0("power_oneside.Trt ",1:K))]
  # plotdata = gather(as.data.frame(unadjustpower), power, estimate, paste0("power_oneside.Trt ",1:K),  factor_key=TRUE)
  # plotdata_power = cbind(plotdata, group = as.numeric(gsub(".+([0-9]+)", "\\1", plotdata$power)))
  #
  # ggplot2::ggplot(data=plotdata_power, aes(x=p_trt, y=estimate, fill=factor(group))) +
  #   geom_bar(data=plotdata_power, aes(x=p_trt, y=estimate, fill=factor(group)),stat="identity",color="black", position=position_dodge())+
  #   geom_text(aes(label=round(estimate,2)), vjust=-0.6, color="black",
  #             position = position_dodge(0.1), size=3.5, angle=30) +
  #   theme_minimal() +
  #   scale_fill_brewer(palette="GnBu", name="Trt") +
  #   ylim(0, 1.15) +
  #   labs(title="Power - One sided proportion test",
  #        subtitle= "Reference to Trt 1 (placebo)",
  #        x ="success rate of trt", y = "power") +
  #   scale_x_continuous(breaks = p_trt)  -> p_oneside
  #
  # ## plot chi-sq test
  #
  # ggplot2::ggplot(data=as.data.frame(outputmat[,c("p_trt","power_chisq")]),
  #        aes(x=p_trt, y=power_chisq)) +
  #   geom_bar(stat="identity",color="black", position=position_dodge(),fill="palegreen4")+
  #   geom_text(aes(label=round(power_chisq,2)), vjust=-0.6, color="black",
  #             position = position_dodge(0.1), size=3.5, angle=30) +
  #   theme_minimal() +
  #   ylim(0, 1.15) +
  #   labs(title="Power - Chi-square test",
  #        x ="success rate of trt", y = "power") +
  #   scale_x_continuous(breaks = p_trt)  -> p_chisq

  ## plot two tests in one plot, x-axis = difference of rates
  powers = as.data.frame(outputmat[,c("p_trt","power_chisq","power_oneside.Trt 2")])
  colnames(powers)[2] <- "Chi-square test"
  colnames(powers)[3] <- "One-sided proportion test"
  plotdata_power = tidyr::gather(powers, test, power,c("Chi-square test","One-sided proportion test"),  factor_key=TRUE)

  ggplot2::ggplot(data= plotdata_power,
                  ggplot2::aes(x=p_trt, y=power, color=test)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() + ggplot2::geom_point() +
    ggplot2::scale_color_brewer(palette="Dark2", name="Tests") +
    ggplot2::ylim(0, 1.) +
    ggplot2::labs(title="Power v.s Treatment response rate",
         subtitle= "Reference to Trt 1 (placebo)",
         x ="success rate of trt", y = "power") +
    ggplot2::theme(legend.position="top",
          legend.justification="right",
          legend.margin=ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10)) +
    ggplot2::scale_x_continuous(breaks = p_trt)  -> combined_power_plot

  result = list(
    Allocation = plotdata_allo,
    Estimation = plotdata_p_est,
    Power = plotdata_power,
    Plot = (p_allocation + p_est) / combined_power_plot
  )
  return(result)
}

