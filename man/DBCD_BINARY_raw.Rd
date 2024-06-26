% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DBCD_BINARY_raw.R
\name{DBCD_BINARY_raw}
\alias{DBCD_BINARY_raw}
\title{Doubly Adaptive Biased Coin Design (Binary Data Frame)}
\usage{
DBCD_BINARY_raw(X.df, group_allo, rho_func_index, rho_func, alpha)
}
\arguments{
\item{X.df}{A data frame of two columns: treatment arm and response value.}

\item{group_allo}{An integer of the size of group allocation. The default is
1.}

\item{rho_func_index}{Supply a number of 1, 2 or 3 indicting the
allocation function to use.
1 = Wei's allocation
2 = Neyman allocation;
3 = Rosenberger allocation.
The default is 3.}

\item{rho_func}{Supply a user-specified allocation function of sample
response rates when rho_func_index is NULL. Default is NULL.}

\item{alpha}{Supply a number indicating the subscripts of the probability
function. The default is 2.}
}
\value{
Code of the arm that the next subject is assigned to.
}
\description{
Allocates patients to one of treatments based on doubly adaptive
biased coin design on individual data.
}
\details{
'DBCD_BINARY_raw' assigns the next subject to a group given the
observed success rates, enrolled subjects and allocation function.
}
\examples{
X.df = data.frame(
ARM = sample(LETTERS[1:3],50,replace = TRUE),
RESPONSE = sample(c(0,1),50,replace = TRUE)
)
DBCD_BINARY_raw(X.df, rho_func_index = 3, alpha=2)

X.df = data.frame(
ARM = sample(LETTERS[1:2],40,replace = TRUE),
RESPONSE = sample(c(0,1),40,replace = TRUE)
)
DBCD_BINARY_raw(
X.df, rho_func_index = NULL,
rho_func = function(x) rev(1-x)/sum(1-x), alpha=2
)

}
