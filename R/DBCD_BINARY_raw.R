#' Doubly Adaptive Biased Coin Design (Binary Data Frame)
#'
#' Allocates patients to one of treatments based on doubly adaptive
#' biased coin design on individual data.
#'
#' @usage DBCD_BINARY_raw(X.df, group_allo, rho_func_index, rho_func, alpha)
#'
#' @param X.df A data frame of two columns: treatment arm and response value.
#' @param group_allo An integer of the size of group allocation. The default is
#' 1.
#' @param rho_func_index Supply a number of 1, 2 or 3 indicting the
#' allocation function to use.
#' 1 = Wei's allocation
#' 2 = Neyman allocation;
#' 3 = Rosenberger allocation.
#' The default is 3.
#' @param rho_func Supply a user-specified allocation function of sample
#' response rates when rho_func_index is NULL. Default is NULL.
#' @param alpha Supply a number indicating the subscripts of the probability
#' function. The default is 2.
#'
#' @importFrom stats na.omit
#'
#' @details 'DBCD_BINARY_raw' assigns the next subject to a group given the
#' observed success rates, enrolled subjects and allocation function.
#'
#' @return Code of the arm that the next subject is assigned to.
#' @export
#'
#' @examples
#' X.df = data.frame(
#' ARM = sample(LETTERS[1:3],50,replace = TRUE),
#' RESPONSE = sample(c(0,1),50,replace = TRUE)
#' )
#' DBCD_BINARY_raw(X.df, rho_func_index = 3, alpha=2)
#'
#' X.df = data.frame(
#' ARM = sample(LETTERS[1:2],40,replace = TRUE),
#' RESPONSE = sample(c(0,1),40,replace = TRUE)
#' )
#' DBCD_BINARY_raw(
#' X.df, rho_func_index = NULL,
#' rho_func = function(x) rev(1-x)/sum(1-x), alpha=2
#' )
#'
DBCD_BINARY_raw = function(X.df,
                           group_allo = 1,
                           rho_func_index = 3,
                           rho_func = NULL,
                           alpha = 2){

  if(!is.data.frame(X.df)) stop("The input must be a data frame with 2 columns: ARM and RESPONSE")
  if(ncol(X.df)==1) stop("The input must be a data frame with 2 columns: ARM and RESPONSE")
  if(ncol(X.df)>2) warning("Only the first two columns are used")
  X.df = na.omit(X.df[,1:2])
  arm = X.df[,1];
  response = X.df[,2]
  if(!is.character(arm) | !length(levels(factor(response)))==2) stop("The input data frame must contain a character variable of ARM and a binary varialbe of RESPONSE")
  if(all(is.null(c(rho_func_index,rho_func)))) stop("Missing allocation function.")

  if(is.null(rho_func_index)) Rho_fun = rho_func else
    Rho_fun = get(paste0("Rho_fun", rho_func_index))

  X.table = table(arm,response)
  N_RK = apply(X.table,1,sum)
  S_RK = X.table[,1]/N_RK
  K = length(S_RK)

  S_RK_adj = (S_RK * N_RK + 1 / K) / (N_RK + 1)
  rho = Rho_fun(S_RK_adj)

  phi = g_fun(N_RK / sum(N_RK), rho, alpha)

  new_assign = sample(rownames(X.table), group_allo, prob = phi,replace = T)

  return(list(new_assign = new_assign, allo_prob = phi))
}
