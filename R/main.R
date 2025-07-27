#' Main MR2G analysis
#'
#' @param b_mat Matrix of SNP-trait effects.
#' @param se_mat Matrix of standard errors corresponding to \code{b_mat}.
#' @param n_vec Vector of per-trait sample sizes.
#' @param rho_mat Trait correlation matrix.
#' @param IV_list List of IVs per trait pair.
#' @param R_list LD / correlation structure list.
#' @param pval_mat Optional matrix of p-values (same shape as \code{b_mat}).
#' @param n_mat Optional matrix of sample sizes per SNP-trait pair (same shape as \code{b_mat}).
#' @param trait_vec Optional vector of trait names.
#' @param IV_screen One of \code{"ScreenAug"}, \code{"ScreenMax"}, or \code{"ScreenPair"}.
#' @param pthres P-value threshold for edge selection in ScreenMax result when \code{IV_screen='ScreenAug'}.
#' @param random_start Number of random starts for optimization in \code{mr_cML_O}.
#' @param num_pert_M Number of perturbations used for ScreenMax/ScreenAug phases.
#' @param num_pert Number of perturbations for final result phase.
#' @param seed Random seed.
#' @param check Logical; Whether or not to exclude perturbation results with spectral radius greater than 1.
#'
#' @return A list (shape depends on IV_screen) of causal graph results.
#'
#' @examples
#' # Load example data
#' data("example", package = "MR2G")
#' # Run MR2G on the example data (with small num_pert for quick execution)
#' result <- MR2G(
#'   b_mat   = b_mat,
#'   se_mat  = se_mat,
#'   n_vec   = n_vec,
#'   rho_mat = rho_mat,
#'   IV_list = IV_list,
#'   R_list  = R_list,
#'   num_pert = 5, num_pert_M = 5  # small test run
#' )
#' str(result)
#'
#'@export
#'
MR2G <- function(b_mat,se_mat,n_vec,rho_mat,IV_list,R_list,
                 pval_mat=NULL,n_mat=NULL,trait_vec=NULL,IV_screen = 'ScreenAug',pthres=0.05/((length(n_vec))^2-length(n_vec)),
                 random_start = 0,num_pert_M = 100,num_pert = 100,seed=0,check=TRUE){
  SP_res = ScreenPair(b_mat = b_mat,se_mat = se_mat,
                      pval_mat = NULL,n_mat = NULL,
                      n_vec = n_vec,IV_list = IV_list,R_list = R_list,rho_mat = rho_mat)
  if(IV_screen %in% c('ScreenAug','ScreenMax')){
    SM_res = ScreenMax(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, pre_screen=SP_res$IJ_snp_list)
    SM_res$DP_mat_list = SP_res$DP_mat_list
    if(IV_screen == 'ScreenMax') {
      message('**** Using ScreenMax for IV selection ****')
      num_pert_M = num_pert
    }else{
      message('**** Using ScreenAug for IV selection ****')
    }
    dp_list_SM = Graph_Perturb_fast(screen_res = SM_res,
                            n_vec = n_vec,
                            rho_mat = rho_mat,
                            random_start = random_start,
                            num_pert = num_pert_M, seed=seed, trait_vec=trait_vec)
    resMax = subset_Graph(dp_list=dp_list_SM, keep_trait='all',check=check)
    if(IV_screen == 'ScreenAug'){
      B = resMax$dp_res$Direct_graph_pval < pthres
      gg = igraph::graph_from_adjacency_matrix(B, mode = 'directed', weighted=FALSE)
      SA_res = ScreenAug(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, pre_screen_pair=SP_res$IJ_snp_list,
                                        pre_screen_max=SM_res$IJ_snp_list, graph = gg)
      SA_res$DP_mat_list = SP_res$DP_mat_list
      dp_list_SA = Graph_Perturb_fast(screen_res = SA_res,
                              n_vec = n_vec,
                              rho_mat = rho_mat,
                              random_start = random_start,
                              num_pert = num_pert_M, seed=seed, trait_vec=trait_vec)
      resAug = subset_Graph(dp_list=dp_list_SA, keep_trait='all',check=check)
      return(list(ScreenAug=lapply(resAug$dp_res, t),ScreenMax=lapply(resMax$dp_res, t)))
    }
    return(lapply(resMax$dp_res, t))
  }else{
    message('**** Using ScreenPair for IV selection ****')
    dp_list_SP = Graph_Perturb_fast(screen_res = SP_res,
                            n_vec = n_vec,
                            rho_mat = rho_mat,
                            random_start = random_start,
                            num_pert = num_pert, seed=seed, trait_vec=trait_vec)
    resPair = subset_Graph(dp_list=dp_list_SP, keep_trait='all',check=check)
    return(lapply(resPair$dp_res, t))
  }

}
