#' @keywords internal
#' @noRd
find_invalid_ivs_for_outcome <- function(graph, exposure, outcome, initial_ivs, Z_connections) {
  # graph: igraph object representing the directed graph (only Y nodes)
  # exposure: The exposure trait
  # outcome: The target trait
  # initial_ivs: A vector of candidate IVs for exp
  # Z_connections: A list of connections from Z to Y nodes (e.g., list(Z1 = c("Y3"), Z2 = c("Y1", "Y4")))

  invalid_ivs <- c()

  for (Z in initial_ivs) {
    # Check if Z has direct connections to any Y node
    if (Z %in% names(Z_connections)) {
      connected_nodes <- Z_connections[[as.character(Z)]]

      for (node in connected_nodes) {
        #        if (node == outcome ||
        #             any(igraph::shortest_paths(graph, from = node, to = outcome)$vpath[[1]] != exposure)) {
        #           invalid_ivs <- c(invalid_ivs, Z)
        #           break
        #       }

        # If Z has a path to the outcome bypassing Y1, mark it as invalid
        g_without_k <- igraph::delete_vertices(graph, exposure)
        # Check if there is a path from i to j in the modified graph
        path_exists <- lengths(suppressWarnings(igraph::shortest_paths(g_without_k, from = node, to = outcome)$vpath)) > 0
        if (node == outcome || path_exists) {
          invalid_ivs <- c(invalid_ivs, Z)
          break
        }
      }
    }
  }

  return(invalid_ivs)
}

#' @keywords internal
#' @noRd
construct_z_connections <- function(effect_matrix, pval_matrix, exposure, sig.cutoff = 5e-8) {
  # effect_matrix: A matrix of marginal effects (rows: IVs, columns: traits)
  # Y1: The trait being considered

  Z_connections <- list()
  effect_matrix = abs(effect_matrix)
  for (iv in rownames(effect_matrix)) {
    # Identify traits where the IV has stronger effects than on Y1
    stronger_index = (effect_matrix[iv, ] > effect_matrix[iv, exposure]) & (pval_matrix[iv,] < sig.cutoff)
    stronger_index[is.na(stronger_index)] = FALSE
    stronger_effects <- colnames(effect_matrix)[stronger_index]

    if (length(stronger_effects) > 0) {
      Z_connections[[iv]] <- stronger_effects
    }
  }

  return(Z_connections)
}

#' @keywords internal
#' @noRd
ScreenPair <- function(b_mat,se_mat,pval_mat,n_mat,n_vec,IV_list,R_list,rho_mat,sig.cutoff=5e-08){
  n_trait = length(n_vec)
  N_combination = n_trait * (n_trait - 1) / 2
  if(length(IV_list)!=N_combination){stop("The length of IV_list must be equal to N_combination!")}
  m_block = length(R_list) ## LD blocks
  if(sum(unlist(lapply(R_list,function(x){nrow(x$R)})))!=nrow(b_mat)){stop("LD matrix must contain all SNPs used in the analysis!")}

  IJ_snp_list = vector("list", N_combination)
  k = 1
  problem_k = NULL
  ### Screening for IVs for each pair of traits ###
  ## i:trait 1; j:trait 2
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      SNP_used = IV_list[[k]]
      b_X = b_mat[,i]
      b_Y = b_mat[,j]
      se_X = se_mat[,i]
      se_Y = se_mat[,j]
      if(is.null(n_mat)){
        n_X = n_vec[i]
        n_Y = n_vec[j]
      }else{
        n_X = n_mat[,i]
        n_Y = n_mat[,j]
      }
      keep_ind = intersect(which(!is.na(b_X)),which(!is.na(b_Y)))
      if(is.null(pval_mat)){
        pvalue.X = pnorm(-abs(b_X/se_X))*2
        pvalue.Y = pnorm(-abs(b_Y/se_Y))*2
      }else{
        pvalue.X = pval_mat[,i]
        pvalue.Y = pval_mat[,j]
      }
      cor_X = b_X / sqrt(b_X^2 + (n_X-2)*se_X^2)
      cor_Y = b_Y / sqrt(b_Y^2 + (n_Y-2)*se_Y^2)

      ind_X = which(pvalue.X<(sig.cutoff))
      ind_Y = which(pvalue.Y<(sig.cutoff))
      # With Screening
      intersect.ind.X.Y = intersect(ind_X,ind_Y)
      ind_X_new = setdiff(ind_X,
                          intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])<
                                              (abs(cor_Y)[intersect.ind.X.Y])])
      ind_Y_new = setdiff(ind_Y,
                          intersect.ind.X.Y[(abs(cor_X)[intersect.ind.X.Y])>
                                              (abs(cor_Y)[intersect.ind.X.Y])])
      ind_X_new = intersect(ind_X_new,keep_ind)
      ind_Y_new = intersect(ind_Y_new,keep_ind)
      ind_X_new_snp = rownames(b_mat)[ind_X_new]
      ind_Y_new_snp = rownames(b_mat)[ind_Y_new]

      ind_X_new1 = ind_X_new[is.element(ind_X_new_snp,SNP_used)]
      ind_Y_new1 = ind_Y_new[is.element(ind_Y_new_snp,SNP_used)]

      IJ_snp_list[[k]]$ind_i_new = ind_X_new1
      IJ_snp_list[[k]]$ind_j_new = ind_Y_new1
      if(length(ind_X_new1)<3){warning(paste0('For trait ', i, ' and trait ',j,', less than 3 IVs were used for trait ', i)); problem_k=c(problem_k,k)}
      if(length(ind_Y_new1)<3){warning(paste0('For trait ', i, ' and trait ',j,', less than 3 IVs were used for trait ', j)); problem_k=c(problem_k,k)}
      k = k + 1
    }
  }
  ### End of screening 1 ###
  ### Screening for LD of b_mat ###
  DP_mat_list = vector("list", m_block)
  for(i in 1:m_block){
    R = R_list[[i]]$R
    match_order = match(rownames(b_mat),rownames(R))
    match_order = match_order[!is.na(match_order)]
    R = R[match_order,match_order,drop=FALSE]
    snp_in_block = R_list[[i]]$snp
    PXR = kronecker(rho_mat,R)
    eigen_decomp = eigen(PXR)
    D = diag(sqrt(zapsmall(eigen_decomp$value,digits=6)))
    DP_mat_list[[i]]$V = eigen_decomp$vector %*% D
    DP_mat_list[[i]]$snp = rownames(R)
  }
  ### End of screening 2 ###

  ### End of screnning 3 ###
  out = list()
  out$IJ_snp_list = IJ_snp_list
  out$DP_mat_list = DP_mat_list
  out$b_mat = b_mat
  out$se_mat = se_mat
  return(out)
}

#' @keywords internal
#' @noRd
ScreenMax <- function(b_mat,se_mat,n_vec,pre_screen, sig.cutoff=5e-08){
  n_trait = length(n_vec)
  N_combination = n_trait * (n_trait - 1) / 2
  a = t(apply(se_mat, 1, function(x){(n_vec-2)*x^2}))
  cor_mat = b_mat / sqrt(b_mat^2 + a)
  rownames(cor_mat) = rownames(b_mat);
  colnames(cor_mat)  = colnames(b_mat);
  pval_mat = pnorm(-abs(b_mat/se_mat))*2
  cor_mat[pval_mat > sig.cutoff] = 0
  if(is.null(rownames(cor_mat))){ rownames(cor_mat) = 1: nrow(cor_mat)};
  if(is.null(colnames(cor_mat))){ colnames(cor_mat) = 1: ncol(cor_mat)};
  IJ_snp_list = vector("list", N_combination)
  IV_max_list = vector("list",n_trait)
  trait_vec = colnames(b_mat)
  if(is.null(trait_vec)){trait_vec = as.character(1:n_trait)}
  if(length(pre_screen)!=N_combination){stop("The length of pre_screen must be equal to N_combination!")}
  for(i in 1:n_trait){
    IV_max_list[[i]] = which(apply(abs(cor_mat),1,which.max)==i)
  }
  k = 1
  ### Screening for IVs for each pair of traits ###
  ## i:trait 1; j:trait 2
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){

      ind_X_new = intersect(pre_screen[[k]]$ind_i_new,IV_max_list[[i]])
      ind_Y_new = intersect(pre_screen[[k]]$ind_j_new,IV_max_list[[j]])
      ind_X_new_snp = rownames(b_mat)[ind_X_new]
      ind_Y_new_snp = rownames(b_mat)[ind_Y_new]

      IJ_snp_list[[k]]$ind_i_new = which(rownames(b_mat) %in% ind_X_new_snp)
      IJ_snp_list[[k]]$ind_j_new = which(rownames(b_mat) %in% ind_Y_new_snp)
      k = k + 1
    }
  }
  ### End of screening 1 ###

  out = list()
  out$IJ_snp_list = IJ_snp_list
  out$b_mat = b_mat
  out$se_mat = se_mat
  return(out)
}

#' @keywords internal
#' @noRd
ScreenAug <- function(b_mat,se_mat,n_vec,pre_screen_pair,pre_screen_max,graph){
  n_trait = length(n_vec)
  N_combination = n_trait * (n_trait - 1) / 2
  a = t(apply(se_mat, 1, function(x){(n_vec-2)*x^2}))
  cor_mat = b_mat / sqrt(b_mat^2 + a)
  rownames(cor_mat) = rownames(b_mat);
  colnames(cor_mat)  = colnames(b_mat);
  pval_mat = pnorm(-abs(b_mat/se_mat))*2
  if(is.null(rownames(cor_mat))){ rownames(cor_mat) = 1: nrow(cor_mat)};
  if(is.null(colnames(cor_mat))){ colnames(cor_mat) = 1: ncol(cor_mat)};
  IJ_snp_list = vector("list", N_combination)
  IV_max_list = vector("list",n_trait)
  trait_vec = colnames(b_mat)
  if(is.null(trait_vec)){trait_vec = as.character(1:n_trait)}
  if(length(pre_screen_pair)!=N_combination){stop("The length of pre_screen must be equal to N_combination!")}
  k = 1
  ### Screening for IVs for each pair of traits ###
  ## i:trait 1; j:trait 2
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){

      ind_X_add = setdiff(pre_screen_pair[[k]]$ind_i_new,pre_screen_max[[k]]$ind_i_new)
      X_invalid_ivs <- find_invalid_ivs_for_outcome(graph = graph, exposure = trait_vec[i], outcome = trait_vec[j],
                                                    initial_ivs = rownames(cor_mat)[ind_X_add],
                                                    Z_connections = construct_z_connections(effect_matrix=cor_mat[ind_X_add,], pval_matrix = pval_mat[ind_X_add,], exposure = trait_vec[i]))
      ind_X_invalid_ivs = which(rownames(b_mat) %in% X_invalid_ivs)
      ind_X_new = union(pre_screen_max[[k]]$ind_i_new, setdiff(ind_X_add,ind_X_invalid_ivs))

      ind_Y_add = setdiff(pre_screen_pair[[k]]$ind_j_new,pre_screen_max[[k]]$ind_j_new)
      Y_invalid_ivs <- find_invalid_ivs_for_outcome(graph = graph, exposure = trait_vec[j], outcome = trait_vec[i],
                                                    initial_ivs = rownames(cor_mat)[ind_Y_add],
                                                    Z_connections = construct_z_connections(cor_mat[ind_Y_add,], pval_matrix = pval_mat[ind_Y_add,],exposure = trait_vec[j]))
      ind_Y_invalid_ivs = which(rownames(b_mat) %in% Y_invalid_ivs)
      ind_Y_new = union(pre_screen_max[[k]]$ind_j_new, setdiff(ind_Y_add,ind_Y_invalid_ivs))

      ind_X_new_snp = rownames(b_mat)[ind_X_new]
      ind_Y_new_snp = rownames(b_mat)[ind_Y_new]

      IJ_snp_list[[k]]$ind_i_new = which(rownames(b_mat) %in% ind_X_new_snp)
      IJ_snp_list[[k]]$ind_j_new = which(rownames(b_mat) %in% ind_Y_new_snp)
      k = k + 1
    }
  }
  ### End of screening 1 ###

  out = list()
  out$IJ_snp_list = IJ_snp_list
  out$b_mat = b_mat
  out$se_mat = se_mat
  return(out)
}

