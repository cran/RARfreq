#' Doubly Adaptive Biased Coin Design (Binary Responses)
#'
#' Allocates patients to one of treatments based on doubly adaptive
#' biased coin design on summary level data.
#'
#' @usage DBCD_BINARY(S_RK, N_RK, group_allo, rho_func_index, rho_func, alpha)
#'
#' @param S_RK A vector of current success rates for each treatment arm.
#' @param N_RK A vector of current enrolled subjects in each arm. The length
#' must be the same as `S_RK`.
#' @param group_allo An integer of the size of group allocation. The default is
#' 1.
#' @param rho_func_index Supply a number of 1, 2 or 3 indicting the
#' allocation function to use.
#' 1 = Wei's allocation;
#' 2 = Neyman allocation;
#' 3 = Rosenberger allocation.
#' The default is 3.
#' @param rho_func Supply a user-specified allocation function of S_RK when
#' rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. The default is 2.
#'
#'
#' @details 'DBCD_BINARY' assigns the next subject to a group given the observed
#' success rates, enrolled subjects and allocation function.
#'
#' @return Number of the arm that the next subject is assigned to.
#' @export
#' @examples
#' # There are 3 arms each of which receives 25 patients.
#' # The observed response rates are 40%, 28% and 60%, respectively.
#' # The following command returns the number of arm that the next patient will
#' # be assigned to.
#' DBCD_BINARY(S_RK = c(0.4, 0.28, 0.6),
#'             N_RK = c(25, 25, 25),
#'             rho_func_index = 3, alpha=2)
#'
#' # Urn allocation
#' DBCD_BINARY(S_RK = c(0.4, 0.3),
#'             N_RK = c(25, 25),
#'             group_allo = 1,
#'             rho_func_index = NULL,
#'             rho_func = function(x) rev(1-x)/sum(1-x),
#'             alpha=2)
#'
DBCD_BINARY = function(S_RK,
                       N_RK,
                       group_allo = 1,
                       rho_func_index = 3,
                       rho_func = NULL,
                       alpha = 2) {
  K = length(S_RK)
  if (length(N_RK) != K)
    stop("S_RK (success rates) and N_RK (sample sizes) should have the same length.")
  if (any(S_RK<0) |any(S_RK>1) )
    stop("S_RK (success rates) should be within 0 and 1.")
  if(all(is.null(c(rho_func_index,rho_func)))) stop("Missing allocation function.")

  if(is.null(rho_func_index)) Rho_fun = rho_func else
    Rho_fun = get(paste0("Rho_fun", rho_func_index))
  S_RK_adj = (S_RK * N_RK + 1 / K) / (N_RK + 1)
  rho = Rho_fun(S_RK_adj)

  phi = g_fun(N_RK / sum(N_RK), rho, alpha) # to be added: set phi=rho

  new_assign = sample(1:K, group_allo, prob = phi,replace = T)

  return(list(new_assign = new_assign, allo_prob = phi))
}


