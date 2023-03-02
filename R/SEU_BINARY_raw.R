#' Sequential Estimation-adjusted Urn Model (Binary Data)
#'
#' Allocates patients to one of treatments based on sequential
#' estimation-adjusted urn model (SEU) on summarized data.
#'
#' @usage SEU_BINARY_raw(x.df, urn_comp, arms, group_allo, add_rule_index,
#' add_rule)
#'
#' @param x.df A data frame of two columns: treatment arm and response value.
#' @param urn_comp A vector of current urn composition.
#' @param arms A vector of arm names. If it is not provided, the arms occurred
#' in x.df will be assumed as all possible arms. Suggest to always assign arms.
#' @param group_allo An integer of the size of group allocation. The default is
#' 1.
#' @param add_rule_index Supply a number of 1, 2 or 3 indicting the
#' addition rules to target allocation functions.
#' 1 = randomized play-the-winner (RPW) rule that targets the urn allocation
#' 2 = the SEU model that targets Neyman allocation;
#' 3 = the SEU model that targets Rosenberger allocation;'
#' 4 = the SEU model that assigns probability of 0.6+1/K to winner at each step.
#' The default is 1.
#' @param add_rule Supply a user-specified addition rules function of x.df and
#' arms when add_rule_index is NULL. Default is NULL.
#'
#' @importFrom stats na.omit
#'
#' @details 'SEU_BINARY_raw' assigns the next subject to a group given the
#' observed data, current urn composition, full list of arm codes,
#' number of group allocation and addition rule function.
#'
#' @return Code of arms that the next group of subjects assigned to and the
#' updated urn composition.
#' @export
#'
#' @examples
#' x.df = data.frame(
#' ARM = sample(LETTERS[1:3],50,replace = TRUE),
#' RESPONSE = sample(c(0,1),50,replace = TRUE)
#' )
#' SEU_BINARY_raw(x.df, urn_comp=c(0,0,0), arms=c("A","B","C"))
#'
#' x.df = data.frame(
#' ARM = sample(LETTERS[1:2],40,replace = TRUE),
#' RESPONSE = sample(c(0,1),40,replace = TRUE)
#' )
#' SEU_BINARY_raw(x.df,
#' urn_comp=c(0,0),
#' arms=c("A","B"),
#' group_allo = 1,
#' add_rule_index = 3)
#'
SEU_BINARY_raw = function(x.df,
                          urn_comp, #Y_n
                          arms=NULL,
                          group_allo = 1,
                          add_rule_index = 1,
                          add_rule = NULL) {
  if(!is.data.frame(x.df)) stop("The input must be a data frame with 2 columns: ARM and RESPONSE")
  if(ncol(x.df)==1) stop("The input must be a data frame with 2 columns: ARM and RESPONSE")
  if(ncol(x.df)>2) warning("Only the first two columns are used")
  x.df = na.omit(x.df[,1:2])
  arm = x.df[,1];
  response = x.df[,2]
  # if(!is.character(arm) | !length(levels(factor(response)))==2) stop("The input data frame must contain a character variable of ARM and a binary variable of RESPONSE")

  if(any(urn_comp<0)) stop("The particles in urn composition cannot be negative.")

  if(all(is.null(c(add_rule_index,add_rule)))) stop("Missing addition rule in SEU.")
  if(is.null(add_rule_index)) Add_rule = add_rule else
    Add_rule = get(paste0("add_rule", add_rule_index))

  if(is.null(arms)) arms = levels(factor(arm))

  urn_comp_update = urn_comp + Add_rule(x.df, arms)
  new_assign = sample(as.character(arms), group_allo, prob = urn_comp_update,replace = T)

  return(list(new_assign = new_assign, urn_comp = urn_comp_update))
}


