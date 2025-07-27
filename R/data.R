#' Example IV list
#'
#' A list of instrumental variables for testing.
#'
#' @format A list with elements:
#' \describe{
#'   \item{trait1}{Character vector of SNPs}
#'   \item{trait2}{Character vector of SNPs}
#' }
#' @examples
#' data("IV_list")
"IV_list"

#' Example R list
#'
#' @format A list containing LD matrices for IVs.
"R_list"

#' Example b_mat matrix
#'
#' @format A matrix of SNP GWAS effects.
"b_mat"

#' Example n_vec
#'
#' @format Numeric vector of sample sizes per trait.
"n_vec"

#' Example rho_mat
#'
#' @format Correlation matrix between GWAS.
"rho_mat"

#' Example se_mat
#'
#' @format Standard error matrix with the same dimensions as b_mat.
"se_mat"
