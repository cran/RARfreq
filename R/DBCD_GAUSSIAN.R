#' Doubly Adaptive Biased Coin Design (Gaussian Responses)
#'
#' Allocates patients to one of treatments based on doubly adaptive
#' biased coin design on summarized data.
#'
#' @usage DBCD_GAUSSIAN(Mean_RK,SD_RK, N_RK, group_allo, rho_func_index,
#' rho_func, alpha)
#'
#' @param Mean_RK A vector of current response means for each treatment arm.
#' @param SD_RK A vector of current response standard deviations for each
#' treatment arm.
#' @param N_RK A vector of current enrolled subjects in each arm. The length
#' must be the same as `Mean_RK`.
#' @param group_allo An integer of the size of group allocation. The default is
#' 1.
#' @param rho_func_index Supply a number of 1 or 2 indicting the
#' allocation function to use.
#' 1 = Zhang-Rosenberger allocation;
#' 2 = Neyman allocation.
#' The default is 2.
#' @param rho_func Supply a user-specified allocation function of Mean_RK and
#' SD_RK when rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. The default is 2.
#'
#'
#' @details 'DBCD_GAUSSIAN' assigns the next subject to a group given the observed
#' success rates, enrolled subjects and allocation function.
#'
#' @return Number of the arm that the next subject is assigned to.
#' @export
#' @examples
#' # There are 2 arms each of which receives 10 patients.
#' # The observed response means are 4.5 and 5, respectively.
#' # The following command returns the number of arm that the next patient will
#' # be assigned to.
#' DBCD_GAUSSIAN(Mean_RK = c(4.5,5),
#' SD_RK = c(1.32,0.72),
#' N_RK = c(10,10),
#' rho_func_index = 2, alpha=2)
#'
#' DBCD_GAUSSIAN(Mean_RK = c(4.5,5),
#' SD_RK = c(1.32,0.72),
#' N_RK = c(10,10),
#' rho_func_index = 1, alpha=2)
#'
DBCD_GAUSSIAN = function(Mean_RK,
                       SD_RK,
                       N_RK,
                       group_allo = 1,
                       rho_func_index = 2,
                       rho_func = NULL,
                       alpha = 2) {
  K = length(Mean_RK)
  if (length(SD_RK) != K | length(N_RK)!=K)
    stop("Mean_RK (mean parameter), SD_RK (standard deviation parameter) and N_RK (sample sizes) should have the same length.")
  if (any(SD_RK<0) )
    stop("SD_RK (standard deviation parameter) cannot be negative.")
  if(all(is.null(c(rho_func_index,rho_func)))) stop("Missing allocation function.")

  if(is.null(rho_func_index)) Rho_fun = rho_func else
    Rho_fun = get(paste0("Rho_fun", rho_func_index,"_Gaussian"))

  rho = Rho_fun(Mean_RK, SD_RK)

  phi = g_fun(N_RK / sum(N_RK), rho, alpha)

  new_assign = sample(1:K, group_allo, prob = phi, replace = T)

  return(list(new_assign = new_assign, allo_prob = phi))
}


