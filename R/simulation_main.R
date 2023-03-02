#' Doubly Adaptive Biased Coin Design with Simulated Data (Binary Responses)
#'
#' Allocates patients to one of treatments based on the doubly adaptive biased
#' coin design with simulated data.
#'
#' @usage
#'
#' simulation_main(n, nstart, p, nstop, replication, group_allo, rho_func_index,
#' rho_func, alpha, sig_level)
#'
#'
#' @param n The number of patients. The default is 500.
#' @param nstart Burn-in sample size of each arm. The default is n/20.
#' @param p A vector containing response probabilities
#' for each treatment arm (where the first element refers to the control arm).
#' The length of p should correspond to the number of arms.
#' The default is p = c(0.3, 0.3, 0.6).
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
#' allocation function to use; default to 3.
#' * 1 = Wei's allocation (Urn allocation);
#' * 2 = Neyman allocation;
#' * 3 (default) = Rosenberger allocation;
#' @param rho_func Supply a user-specified allocation function of S_RK when
#' rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. The default is 2.
#' @param sig_level Significant level (one-sided). The default is 0.05.
#'
#'
#' @export
#'
#' @details 'simulation_main' can sample response and adaptively randomize
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
#' simulation_main(n = 500,
#' nstart = round(500 / 20),
#' p = c(0.3, 0.3, 0.6),
#' nstop=c(500,500,500),
#' replication = 5,
#' group_allo = 1,
#' rho_func_index = 3,
#' alpha = 2,
#' sig_level = 0.05
#' )
#'

simulation_main = function(n = 500,
                           nstart = NULL,
                           p = c(0.3, 0.3, 0.6),
                           nstop = NULL,
                           replication = 100,
                           group_allo = 1,
                           rho_func_index = 3,
                           rho_func = NULL,
                           alpha = 2,
                           sig_level = 0.05) {

  if(all(is.null(c(rho_func_index,rho_func)))) stop("Missing allocation function.")
  if(is.null(rho_func_index)) Rho_fun = rho_func else
    Rho_fun = get(paste0("Rho_fun", rho_func_index))

  K = length(p)
  if(is.null(nstart)) nstart = round(n/20)
  nstart = round(nstart)
  if(K*nstart>n) stop("No subject to allocate.")
  if(is.null(nstop)) nstop = rep(n, K) else if(length(nstop)!=K)
    stop("The length stopping cap should be the same as number of groups." )

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
    X_initial = sapply(p, function(x)
      stats::rbinom(nstart, 1, x))
    S_RK = apply(X_initial, 2, function(x) mean(c(1/K,x),na.rm=T)) #adjust
    N_RK = rep(nstart, K)

    new_assign = rep(0, (n - K * nstart))
    new_obs = rep(0, (n - K * nstart))

    ###----- assign subjects and draw observations -----###
    X_new = X_initial
    for(j in seq(length(seq_allo))[-1]){
      mydbcd = DBCD_BINARY(S_RK, N_RK,group_allo = seq_allo[j]-seq_allo[j-1],
                           rho_func_index, rho_func, alpha)
      mynext = mydbcd$new_assign
      myobs = sapply(p[mynext],function(x) stats::rbinom(1, 1, x))
      new_assign[(seq_allo[j-1]+1):seq_allo[j] - nstart*K] = mynext
      new_obs[(seq_allo[j-1]+1):seq_allo[j] - nstart*K] = myobs

      X_tmp = matrix(NA, seq_allo[j]-seq_allo[j-1], K)
      X_tmp[cbind(seq(seq_allo[j]-seq_allo[j-1]),mynext)] = myobs
      X_new = rbind(X_new,X_tmp)

      N_RK = apply(X_new,2,function(x) length(na.omit(x)))
      S_RK = apply(X_new,2,function(x) mean(c(1/K,x),na.rm=T)) #adjust

      if(any(N_RK>=nstop)) break

    }
    mydata.df = data.frame(arm = c(rep(1:K,each=nstart),new_assign), measure = c(c(X_initial),new_obs))


    N_Rk_Rrepli[m, ] = N_RK
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
