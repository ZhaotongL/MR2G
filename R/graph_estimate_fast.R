#' Precompute per-pair meta-info
#' @keywords internal
#' @noRd
make_pair_info <- function(n_trait, IJ_snp_list, rho_mat, n_vec) {
  pairs <- utils::combn(n_trait, 2, simplify = FALSE)
  lapply(seq_along(pairs), function(k) {
    i <- pairs[[k]][1]; j <- pairs[[k]][2]
    list(
      i = i, j = j,
      ind_i_new = IJ_snp_list[[k]]$ind_i_new,
      ind_j_new = IJ_snp_list[[k]]$ind_j_new,
      rho = rho_mat[i, j],
      n = min(n_vec[i], n_vec[j])
    )
  })
}

#' Fast Graph Estimation (single perturbation)
#'
#' @param b_mat perturbed effect matrix
#' @param se_mat standard error matrix
#' @param pair_info list from \code{make_pair_info}
#' @param random_start passed to mr_cML_O
#' @param seed random seed
#'
#' @return list: obs_graph, obs_graph_se, obs_graph_pval
#' @export
Graph_Estimate_fast <- function(b_mat, se_mat, pair_info,
                                random_start = 10, seed = 0) {

  if (seed != 0) set.seed(seed)

  n_trait <- ncol(b_mat)
  obs_graph      <- matrix(1, nrow = n_trait, ncol = n_trait)
  obs_graph_se   <- matrix(0, nrow = n_trait, ncol = n_trait)
  obs_graph_pval <- matrix(0, nrow = n_trait, ncol = n_trait)

  for (info in pair_info) {
    i <- info$i; j <- info$j

    ItoJ <- mr_cML_O(
      b_exp = b_mat[info$ind_i_new, i],
      b_out = b_mat[info$ind_i_new, j],
      se_exp = se_mat[info$ind_i_new, i],
      se_out = se_mat[info$ind_i_new, j],
      n = info$n,
      rho = info$rho,
      random_start = random_start,
      random_seed = seed
    )
    JtoI <- mr_cML_O(
      b_exp = b_mat[info$ind_j_new, j],
      b_out = b_mat[info$ind_j_new, i],
      se_exp = se_mat[info$ind_j_new, j],
      se_out = se_mat[info$ind_j_new, i],
      n = info$n,
      rho = info$rho,
      random_start = random_start,
      random_seed = seed
    )

    obs_graph[i, j]      <- ItoJ$BIC_theta
    obs_graph[j, i]      <- JtoI$BIC_theta
    obs_graph_se[i, j]   <- ItoJ$BIC_se
    obs_graph_se[j, i]   <- JtoI$BIC_se
    obs_graph_pval[i, j] <- ItoJ$BIC_p
    obs_graph_pval[j, i] <- JtoI$BIC_p
  }

  list(
    obs_graph      = obs_graph,
    obs_graph_se   = obs_graph_se,
    obs_graph_pval = obs_graph_pval
  )
}
