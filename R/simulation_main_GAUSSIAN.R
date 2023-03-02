#' Doubly Adaptive Biased Coin Design with Simulated Data (Gaussian Responses)
#'
#' Allocates patients to one of treatments based on the doubly adaptive biased
#' coin design with simulated data.
#'
#' @usage
#' simulation_main_GAUSSIAN(n, nstart, mu, sd, nstop, replication, group_allo,
#'  rho_func_index, rho_func, alpha, sig_level)
#'
#'
#' @param n The number of patients. The default is 500.
#' @param nstart Burn-in sample size of each arm. The default is n/20.
#' @param mu A vector of mean response for each treatment arm
#' (where the first element refers to the control arm).
#' The length of mu should correspond to the number of arms.
#' The default is mu = c(4.5,5).
#' @param sd A vector of response standard deviations for each treatment arm.
#' (where the first element refers to the control arm).
#' The length of sd should correspond to the number of arms.
#' The default is sd = c(1.32, 0.72).
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
#' function. The default is 2.
#' @param sig_level Significant level. The default is 0.05.
#'
#'
#' @export
#'
#' @details 'simulation_main_GAUSSIAN' can sample response and adaptively
#' randomize subjects group by group.
#'
#' @return
#' \itemize{
#'   \item allocation_mean - Average of allocation in each arm based on `replication` repeats.
#'   \item allocation_sd - Standard deviation of allocation in each arm based on `replication` repeats.
#'   \item SS_mean - Average of sample size in each arm based on `replication` repeats.
#'   \item SS_sd - Standard deviation of sample size in each arm based on `replication` repeats.
#'   \item power_aov - Average power of ANOVA test.
#'   \item power_oneside - Average power for each of the k-th arm to perform one-sided Welch T-test against H0: mu_1>mu_k without multiplicity adjustment.
#'   \item mu_estimate_mean - Average of estimated response mean `mu`.
#'   \item sd_estimate_mean - Average of estimated response standard deviation `sd`.
#'   \item mu_estimate_sd - Standard deviation of estimated response mean `mu`.
#'   \item sd_estimate_sd - Standard deviation of estimated response standard deviation `sd`.
#' }
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
#' ## Default method
#' simulation_main_GAUSSIAN(
#' n = 500,
#' nstart = round(500 / 20),
#' mu = c(4.5,5),
#' sd = c(1.32,0.72),
#' nstop=c(500,500),
#' replication = 5,
#' group_allo = 1,
#' rho_func_index = 2,
#' rho_func = NULL,
#' alpha = 2,
#' sig_level = 0.05
#' )
#'
simulation_main_GAUSSIAN = function(n = 500,
                                    nstart = NULL,
                                    mu = c(4.5,5),
                                    sd = c(1.32,0.72),
                                    nstop = NULL,
                                    replication = 100,
                                    group_allo = 1,
                                    rho_func_index = 2,
                                    rho_func = NULL,
                                    alpha = 2,
                                    sig_level = 0.05) {

  K = length(mu)
  if (length(sd) != K )
    stop("mu (mean parameter) and sd (standard deviation parameter) should have the same length.")
  if(is.null(nstart)) nstart = round(n/20)
  nstart = round(nstart)
  if(K*nstart>=n) stop("No subject to allocate.")
  if (any(sd<0) )
    stop("sd (standard deviation parameter) cannot be negative.")
  if(is.null(nstop)) nstop = rep(n, K) else if(length(nstop)!=K)
    stop("The length stopping cap should be the same as number of groups." )

  if(any(!is.wholenumber(group_allo))) stop("The group size of allocation must be integer.")
  if(length(group_allo) == 1) seq_allo = seq(from = nstart*K, to = n, by = group_allo)[-1] else
    seq_allo = group_allo
  if(any(seq_allo <= nstart*K)) warning("The allocation before burn-in stage will be skipped.")
  seq_allo = seq_allo[seq_allo >= nstart*K & seq_allo<n]
  seq_allo = sort(unique(c(seq_allo,nstart*K,n)))


  N_Rk_Rrepli = matrix(0, replication, K) #recording total assignment of each group
  Mean_Rk_Rrepli = matrix(0, replication, K)
  SD_Rk_Rrepli = matrix(0, replication, K)

  pval_aov = rep(1, replication)
  pval_oneside = matrix(1, replication, K)
  colnames(pval_oneside) = paste0("Trt ", 1:K)

  for (m in 1:replication) {
    X_initial = mapply(function(x,y) stats::rnorm(nstart, mean=x, sd=y), mu, sd)
    X_new = X_initial
    Mean_RK = apply(X_initial, 2, mean)
    SD_RK = apply(X_initial, 2, sd)
    N_RK = rep(nstart, K)

    new_assign = rep(0, (n - K * nstart))
    new_obs = rep(0, (n - K * nstart))

    for (j in seq(length(seq_allo))[-1]) {
      mydbcd = DBCD_GAUSSIAN(Mean_RK = Mean_RK, SD_RK = SD_RK, N_RK = N_RK,
                             group_allo = seq_allo[j]-seq_allo[j-1],
                             rho_func_index = rho_func_index,
                             rho_func = rho_func,
                             alpha = alpha)
      mynext = mydbcd$new_assign
      myobs = apply(cbind(mu[mynext],sd[mynext]),1, function(x) stats::rnorm(1, x[1], x[2]))
      new_assign[(seq_allo[j-1]+1):seq_allo[j] - nstart*K] = mynext
      new_obs[(seq_allo[j-1]+1):seq_allo[j] - nstart*K] = myobs

      X_tmp = matrix(NA, seq_allo[j]-seq_allo[j-1], K)
      X_tmp[cbind(seq(seq_allo[j]-seq_allo[j-1]),mynext)] = myobs
      X_new = rbind(X_new,X_tmp)

      Mean_RK = apply(X_new, 2, function(x) mean(x,na.rm=T))
      SD_RK = apply(X_new, 2, function(x) sd(x,na.rm=T))
      N_RK = apply(X_new,2,function(x) length(na.omit(x)))

      if(any(N_RK>=nstop)) break
    }

    N_Rk_Rrepli[m, ] = N_RK
    Mean_Rk_Rrepli[m, ] = Mean_RK
    SD_Rk_Rrepli[m, ] = SD_RK

    pval_aov[m] = test_aov(X_new,sig_level)$p.value
    pval_oneside[m, 2:K] = test_WelchT(X_new,sig_level,var.equal=F)$p.value

  }

  allocation_mean = base::apply(N_Rk_Rrepli / n, 2, base::mean)
  allocation_sd = base::apply(N_Rk_Rrepli / n, 2, stats::sd)
  SS_mean = base::apply(N_Rk_Rrepli, 2, base::mean)
  SS_sd = base::apply(N_Rk_Rrepli, 2, stats::sd)
  power_aov = base::mean(pval_aov <= sig_level)
  power_oneside = base::apply(pval_oneside, 2, function(x)
    base::mean(x <= sig_level))
  mu_estimate_mean = base::apply(Mean_Rk_Rrepli, 2, base::mean)
  sd_estimate_mean = base::apply(SD_Rk_Rrepli, 2, base::mean)
  mu_estimate_sd = base::apply(Mean_Rk_Rrepli, 2, sd)
  sd_estimate_sd = base::apply(SD_Rk_Rrepli, 2, sd)

  rval = list(
    allocation_mean = allocation_mean,
    allocation_sd = allocation_sd,
    SS_mean = SS_mean,
    SS_sd = SS_sd,
    power_aov = power_aov,
    power_oneside = power_oneside,
    mu_estimate_mean = mu_estimate_mean,
    sd_estimate_mean = sd_estimate_mean,
    mu_estimate_sd = mu_estimate_sd,
    sd_estimate_sd = sd_estimate_sd
  )

  return(rval)
}
