#' Doubly Adaptive Biased Coin Design (Gaussian Responses)
#'
#' Allocates patients to one of treatments based on doubly adaptive
#' biased coin design on summarized data.
#'
#' @usage DBCD_GAUSSIAN_raw(X.df, group_allo, rho_func_index, rho_func, alpha)
#'
#' @param X.df A data frame of two columns: treatment arm and response value.
#' treatment arm.
#' @param group_allo An integer of the size of group allocation. The default is
#' 1.
#' @param rho_func_index Supply a number of 1, 2 or 3 indicting the
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
#' @importFrom magrittr "%>%"
#' @importFrom stats na.omit
#' @import dplyr
#'
#' @details 'DBCD_GAUSSIAN' assigns the next subject to a group given the observed
#' success rates, enrolled subjects and allocation function.
#'
#' @return Number of the arm that the next subject is assigned to.
#' @export
#' @examples
#' X.df = data.frame(
#' ARM = sample(LETTERS[1:2],50,replace = TRUE),
#' RESPONSE = rnorm(50)
#' )
#'
#' DBCD_GAUSSIAN_raw(X.df)
#'
DBCD_GAUSSIAN_raw = function(X.df,
                             group_allo = 1,
                             rho_func_index = 2,
                             rho_func = NULL,
                             alpha = 2) {
  if(!is.data.frame(X.df)) stop("The input must be a data frame with 2 columns: ARM and RESPONSE")
  if(ncol(X.df)==1) stop("The input must be a data frame with 2 columns: ARM and RESPONSE")
  if(ncol(X.df)>2) warning("Only the first two columns are used")
  X.df = na.omit(X.df[,1:2])
  arm = X.df[,1];
  response = X.df[,2]
  if(!is.character(arm) | !is.numeric(response)) stop("The input data frame must contain a character variable of ARM and a numeric varialbe of RESPONSE")
  if(all(is.null(c(rho_func_index,rho_func)))) stop("Missing allocation function.")

  X.summary = data.frame(arm,response) %>%
    dplyr::group_by(arm) %>%
    dplyr::summarize(Mean_RK = mean(response), SD_RK=sd(response), N_RK=dplyr::n())
  Mean_RK= X.summary$Mean_RK
  SD_RK = X.summary$SD_RK
  N_RK = X.summary$N_RK
  K = length(Mean_RK)

  if(is.null(rho_func_index)) Rho_fun = rho_func else
    Rho_fun = get(paste0("Rho_fun", rho_func_index,"_Gaussian"))

  rho = Rho_fun(Mean_RK, SD_RK)

  phi = g_fun(N_RK / sum(N_RK), rho, alpha)

  new_assign = sample(as.character(unique(arm)), group_allo, prob = phi,replace = T)
  return(list(new_assign = new_assign, allo_prob = phi))
}


