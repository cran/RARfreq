#' Sequential Estimation-adjusted Urn Model with Simulated Data (Binary Data)
#'
#' Allocates patients to one of treatments based on sequential
#' estimation-adjusted urn model (SEU) with simulated data.
#'
#' @usage
#' SEU_simulation_main(n, nstart, p, urn_comp, nstop, replication, group_allo,
#'  add_rule_index, add_rule, add_rule_full, sig_level)
#'
#'
#' @param n The number of patients. The default is 500.
#' @param nstart Burn-in sample size of each arm. The default is n/20.
#' @param p A vector containing response probabilities
#' for each treatment arm (where the first element refers to the control arm).
#' The length of p should correspond to the number of arms.
#' The default is p = c(0.3, 0.3, 0.6).
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
#' @export
#'
#' @details 'SEU_simulation_main' can sample response and adaptively randomize
#' subjects group by group.
#'
#' @return
#' \itemize{
#'   \item allocation_mean - Average of allocation in each arm based on `replication` repeats
#'   \item allocation_sd - Standard deviation of allocation in each arm based on `replication` repeats
#'   \item SS_mean - Average of sample size in each arm based on `replication` repeats
#'   \item SS_sd - Standard deviation of sample size in each arm based on `replication` repeats
#'   \item power_chisq - Average power of chi-square test.
#'   \item power_oneside - Average power for each of the k-th arm to perform one-sided test against H0: p_1>p_k without multiplicity adjustment
#'   \item p_estimate_mean - Average of estimated success rate p
#'   \item p_estimate_sd - Standard deviation of estimated success rate p
#' }
#'
#' @importFrom stats sd
#' @import dplyr
#' @rawNamespace import(data.table, except = c(last, first, between))
#' @examples
#'
#' ## Default method
#' SEU_simulation_main(n = 500,
#' nstart = round(500 / 20),
#' p = c(0.3, 0.3, 0.6),
#' nstop=c(500,500,500),
#' urn_comp = c(0,0,0),
#' replication = 5,
#' group_allo = 1,
#' add_rule_index = 1,
#' add_rule_full = FALSE,
#' sig_level = 0.05
#' )
#'

SEU_simulation_main = function(n = 500,
                               nstart = NULL,
                               p = c(0.3, 0.3, 0.6),
                               urn_comp = NULL,
                               nstop = NULL,
                               replication = 100,
                               group_allo = 1,
                               add_rule_index = 1,
                               add_rule = NULL,
                               add_rule_full = NULL,
                               sig_level = 0.05) {
  ARM <- RESPONSE <- NULL
  # if(all(is.null(c(add_rule_index,add_rule)))) stop("Missing addition rule in SEU.")
  # if(is.null(add_rule_index)) add_rule = add_rule else
  #   Add_rule = get(paste0("add_rule", add_rule_index))
  if(is.null(add_rule_full)){
    add_rule_full = TRUE
    if(!is.null(add_rule_index)) if(add_rule_index==1) add_rule_full = FALSE
  }

  K = length(p)
  if(is.null(nstart)) nstart = round(n/20)
  nstart = round(nstart)
  if(K*nstart>n) stop("No subject to allocate.")
  if(is.null(nstop)) nstop = rep(n, K) else if(length(nstop)!=K)
    stop("The length stopping cap should be the same as number of groups." )
  if(is.null(urn_comp)) urn_comp=rep(0,K) else if(length(urn_comp)!=K)
    stop("The length urn composition should be the same as number of groups." )

  if(any(!is.wholenumber(group_allo))) stop("The group size of allocation must be integer.")
  if(length(group_allo) == 1) seq_allo = seq(from = nstart*K, to = n, by = group_allo)[-1] else
    seq_allo = group_allo
  if(any(seq_allo <= nstart*K)) warning("The allocation before burn-in stage will be skipped.")
  seq_allo = seq_allo[seq_allo >= nstart*K & seq_allo<n]
  seq_allo = sort(unique(c(seq_allo,nstart*K,n)))


  N_Rk_Rrepli = matrix(0, replication, K) #recording total assignment of each group
  S_Rk_Rrepli = matrix(0, replication, K)

  pval_chisq = rep(1, replication)
  pval_oneside = matrix(1, replication, K)
  colnames(pval_oneside) = paste0("Trt ", 1:K)

  for (m in 1:replication) {
    X_initial = sapply(p, function(x) stats::rbinom(nstart, 1, x))
    X.df = data.frame(
     ARM = rep(1:K, each=nstart),
     RESPONSE = as.vector(X_initial)
     )
    # S_RK = apply(X_initial, 2, mean)
    # N_RK = rep(nstart, K)
    new_assign = rep(0, (n - K * nstart))
    new_obs = rep(0, (n - K * nstart))

    ###----- assign subjects and draw observations -----###
    Urn_comp = urn_comp
    X_new = X.df
    X_plug = X.df
    for(j in seq(length(seq_allo))[-1]){
      myseu = SEU_BINARY_raw(X_plug, Urn_comp, arms = 1:K,
                             group_allo = seq_allo[j]-seq_allo[j-1],
                             add_rule_index = add_rule_index, add_rule= add_rule)
      mynext = as.numeric(myseu$new_assign)
      myobs = sapply(p[mynext],function(x) stats::rbinom(1, 1, x))
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
    myseu = SEU_BINARY_raw(X_plug, Urn_comp, arms = 1:K,
                           group_allo = 1,
                           add_rule_index = add_rule_index, add_rule= add_rule)
    Urn_comp = myseu$urn_comp

    # mydata.df = data.frame(arm = c(rep(1:K,each=nstart),new_assign), measure = c(c(X_initial),new_obs))
    X_summary = X_new %>% group_by(ARM) %>%
      dplyr::summarise(mean.adj = base::mean(c(RESPONSE,1/K)), n=n())
    S_RK = X_summary$mean.adj

    N_Rk_Rrepli[m, ] = N_RK #same as X_summary$n
    S_Rk_Rrepli[m, ] = S_RK

    pval_chisq[m] = test_chisq(S_RK, N_RK, sig_level)$p.value
    pval_oneside[m, 2:K] = test_oneside(S_RK, N_RK, sig_level)$p.value

  }

  allocation_mean = apply(N_Rk_Rrepli / rowSums(N_Rk_Rrepli), 2, base::mean)
  allocation_sd = apply(N_Rk_Rrepli / rowSums(N_Rk_Rrepli), 2, stats::sd)
  SS_mean = apply(N_Rk_Rrepli, 2, base::mean)
  SS_sd = apply(N_Rk_Rrepli, 2, stats::sd)
  power_chisq = mean(pval_chisq <= sig_level)
  power_oneside = apply(pval_oneside, 2, function(x)
    mean(x <= sig_level))
  p_estimate_mean = apply(S_Rk_Rrepli, 2, mean)
  p_estimate_sd = apply(S_Rk_Rrepli, 2, sd)

  rval = list(
    allocation_mean = allocation_mean,
    allocation_sd = allocation_sd,
    SS_mean = SS_mean,
    SS_sd = SS_sd,
    power_chisq = power_chisq,
    power_oneside = power_oneside,
    p_estimate_mean = p_estimate_mean,
    p_estimate_sd = p_estimate_sd
  )

  return(rval)
}
