#' Data perturbation of GWAS summary statistics
#' @keywords internal
#' @noRd
Generate_Perturb <- function(b_mat,se_mat,rho_mat,DP_mat_list=NULL,seed=0){
  if(seed!=0){set.seed(seed)}
  n_trait = ncol(b_mat)
  m_used = nrow(b_mat)
  e_mat_dp = matrix(0,ncol=n_trait,nrow=m_used)
  if(!is.null(DP_mat_list)){
    for(i in 1:length(DP_mat_list)){
      V = DP_mat_list[[i]]$V
      X = rnorm(nrow(V))
      e_vec = V %*% X
      e_mat = matrix(e_vec,ncol=n_trait,byrow=FALSE)
      snp_ind = which(is.element(rownames(b_mat),DP_mat_list[[i]]$snp))
      e_mat_dp[snp_ind,] = se_mat[snp_ind,,drop=FALSE] * e_mat
    }
  }else{
    e_mat = MASS::mvrnorm(n=m_used, mu=rep(0,n_trait), Sigma=rho_mat);
    e_mat_dp = se_mat * e_mat;
  }
  b_mat_dp = b_mat + e_mat_dp

  return(b_mat_dp)
}

#' Parallel graph perturbation (fast)
#'
#' @param n_vec See \code{MR2G()}.
#' @param rho_mat See \code{MR2G()}.
#' @param screen_res A list that must contain \code{b_mat}, \code{se_mat},
#'        \code{IJ_snp_list}, and optionally \code{DP_mat_list}.
#' @param num_pert Number of perturbations.
#' @param random_start See \code{MR2G()}.
#' @param seed Random seed.
#' @param trait_vec Optional trait names.
#'
#' @return List containing lists of \code{obs_graph}, \code{obs_graph_se}, \code{obs_graph_pval}, and \code{trait_vec}.
#' @export
Graph_Perturb_fast <- function(n_vec, rho_mat, screen_res,
                               num_pert = 100,
                               random_start = 10,
                               seed = 0,
                               trait_vec = NULL) {

  if (is.null(trait_vec)) {
    trait_vec <- colnames(screen_res$b_mat)
  }

  set.seed(seed)

  n_trait   <- length(n_vec)
  pair_info <- make_pair_info(n_trait, screen_res$IJ_snp_list, rho_mat, n_vec)

  dp_blocks <- NULL
  if (!is.null(screen_res$DP_mat_list)) {
    rn <- rownames(screen_res$b_mat)
    dp_blocks <- lapply(screen_res$DP_mat_list, function(d) {
      list(V = d$V, snp_ind = match(d$snp, rn))
    })
  }

  workers    <- future::nbrOfWorkers()
  chunk_size <- ceiling(num_pert / max(workers, 1))

  progressr::with_progress({
    p <- progressr::progressor(steps = num_pert)

    out_list <- future.apply::future_lapply(
      X = seq_len(num_pert),
      future.seed = TRUE,
      future.chunk.size = chunk_size,
      FUN = function(it) {
        b_mat_dp <- Generate_Perturb(
          b_mat     = screen_res$b_mat,
          se_mat    = screen_res$se_mat,
          rho_mat   = rho_mat,
          DP_mat_list = screen_res$DP_mat_list
        )
        res <- Graph_Estimate_fast(
          b_mat        = b_mat_dp,
          se_mat       = screen_res$se_mat,
          pair_info    = pair_info,
          random_start = random_start
        )
        p()
        res
      }
    )
  })

  list(
    obs_graph_list      = lapply(out_list, `[[`, "obs_graph"),
    obs_graph_se_list   = lapply(out_list, `[[`, "obs_graph_se"),
    obs_graph_pval_list = lapply(out_list, `[[`, "obs_graph_pval"),
    trait_vec           = trait_vec
  )
}
