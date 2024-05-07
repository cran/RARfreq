#' Sequential Estimation-adjusted Urn Model with Simulated Data (Gaussian
#' Responses)
#'
#' Allocates patients to one of treatments based on sequential
#' estimation-adjusted urn model (SEU) with simulated data.
#'
#' @usage
#' SEU_simulation_main_GAUSSIAN(n, nstart, mu, sd, urn_comp, nstop,
#'  replication, group_allo, add_rule_index, add_rule,
#'  add_rule_full, sig_level)
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
#' @param urn_comp A vector of current urn composition. The default is NULL,
#' which indicates no ball in the urn.
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
#' @param add_rule_index Supply a number of 1 or 2 indicting the
#' addition rules to target allocation functions.
#' 1 = the SEU model targeting Neyman allocation;
#' 2 = the SEU model that assigns probability of 0.6+1/K to winner at each step.
#' The default is 1.
#' @param add_rule Supply a user-specified addition rules function of x.df and
#' arms when add_rule_index is NULL. Default is NULL. (See SEU_GAUSSIAN_raw for
#' details on x.df and arms.)
#' @param add_rule_full Indicator of reference data for updating addition rule.
#' If TRUE, the addition rule is updated by full observation at each group
#' allocation. If FALSE,the addition rule is updated by each group observation.
#' The default is TRUE.
#' @param sig_level Significant level (one-sided). The default is 0.05.
#'
#'
#' @export
#'
#' @details 'SEU_simulation_main_GAUSSIAN' can sample response and adaptively
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
#' @importFrom stats model.matrix
#' @import patchwork
#' @import dplyr
#' @import ggplot2
#'
#' @examples
#'
#' ## Default method
#' SEU_simulation_main_GAUSSIAN(
#' n = 50,
#' nstart = round(50 / 20),
#' mu = c(4.5,5),
#' sd = c(1.32,0.72),
#' urn_comp = c(0,0),
#' nstop=c(50,50),
#' replication = 5,
#' group_allo = 1,
#' add_rule_index = 1,
#' add_rule = NULL,
#' add_rule_full = TRUE,
#' sig_level = 0.05
#' )
#'
SEU_simulation_main_GAUSSIAN = function(n = 500,
                                        nstart = NULL,
                                        mu = c(4.5,5),
                                        sd = c(1.32,0.72),
                                        urn_comp = NULL,
                                        nstop = NULL,
                                        replication = 100,
                                        group_allo = 1,
                                        add_rule_index = 1,
                                        add_rule = NULL,
                                        add_rule_full = NULL,
                                        sig_level = 0.05) {
  ARM <- RESPONSE <- NULL
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
  if(is.null(urn_comp)) urn_comp=rep(0,K) else if(length(urn_comp)!=K)
    stop("The length urn composition should be the same as number of groups." )
  if(is.null(add_rule_full))  add_rule_full = TRUE


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
    X.df = data.frame(
      ARM = rep(1:K, each=nstart),
      RESPONSE = as.vector(X_initial)
    )

    new_assign = rep(0, (n - K * nstart))
    new_obs = rep(0, (n - K * nstart))
    ###----- assign subjects and draw observations -----###
    Urn_comp = urn_comp
    X_new = X.df
    X_plug = X.df
    for (j in seq(length(seq_allo))[-1]) {
      myseu = SEU_GAUSSIAN_raw(x.df = X_plug,
                               urn_comp = Urn_comp,
                               arms = 1:K,
                               group_allo = seq_allo[j]-seq_allo[j-1],
                               add_rule_index = add_rule_index, add_rule= add_rule)
      mynext = as.numeric(myseu$new_assign)
      myobs = apply(cbind(mu[mynext],sd[mynext]),1, function(x) stats::rnorm(1, x[1], x[2]))
      new_assign[(seq_allo[j-1]+1):seq_allo[j] - nstart*K] = mynext
      new_obs[(seq_allo[j-1]+1):seq_allo[j] - nstart*K] = myobs
      X_tmp = data.frame(ARM = mynext, RESPONSE = myobs)
      X_new = rbind(X_new,X_tmp)
      if(add_rule_full) X_plug = X_new else X_plug = X_tmp
      Urn_comp = myseu$urn_comp
      N_RK = table(X_new$ARM)

      if(any(N_RK>=nstop)) break
    }

    #final update for the urn
    myseu = SEU_GAUSSIAN_raw(X_plug, Urn_comp, arms = 1:K,
                             group_allo = 1,
                             add_rule_index = add_rule_index, add_rule= add_rule)
    Urn_comp = myseu$urn_comp

    X_summary = X_new %>% group_by(ARM) %>%
      dplyr::summarise(mean = base::mean(RESPONSE), sd = stats::sd(RESPONSE), n=n())
    Mean_RK = X_summary$mean
    SD_RK = X_summary$sd

    N_Rk_Rrepli[m, ] = N_RK
    Mean_Rk_Rrepli[m, ] = Mean_RK
    SD_Rk_Rrepli[m, ] = SD_RK

    pval_aov[m] = stats::anova(lm(RESPONSE~ARM, data=X_new))$Pr[1]
    pval_oneside[m, 2:K] = test_WelchT_long(X_new,sig_level,var.equal=F)$p.value

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
  mu_estimate_sd = base::apply(Mean_Rk_Rrepli, 2, stats::sd)
  sd_estimate_sd = base::apply(SD_Rk_Rrepli, 2, stats::sd)

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
