#' @keywords internal
#' @noRd
Alg1 <- function(M_mat,thres = 1e-8){
  p = nrow(M_mat)
  G_tot = M_mat - diag(p)
  G_dir = matrix(0, nrow=p, ncol = p);
  a = diag(solve(M_mat))
  W = diag(a) %*% M_mat
  G_dir = diag(p) - solve(W)

  G_dir[abs(G_dir)<thres] = 0
  return(G_dir)
}

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

#' @keywords internal
#' @noRd
Graph_Estimate <- function(b_mat,se_mat,n_vec,rho_mat,IJ_snp_list,random_start=10,seed=0){
  n_trait = length(n_vec)
  obs_graph = matrix(1,nrow=n_trait,ncol=n_trait)
  obs_graph_pval = obs_graph_se = matrix(0,nrow=n_trait,ncol=n_trait)
  k = 1
  for(i in 1:(n_trait-1)){
    for(j in (i+1):n_trait){
      ind_i_new = IJ_snp_list[[k]]$ind_i_new
      ind_j_new = IJ_snp_list[[k]]$ind_j_new
      rho_ij = rho_mat[i,j]
      ItoJ_cML_O_res = mr_cML_O(b_exp=b_mat[ind_i_new,i],
                                b_out=b_mat[ind_i_new,j],
                                se_exp=se_mat[ind_i_new,i],
                                se_out=se_mat[ind_i_new,j],
                                n = min(n_vec[i],n_vec[j]),
                                rho = rho_ij,
                                random_start = random_start,random_seed=seed)
      JtoI_cML_O_res = mr_cML_O(b_exp=b_mat[ind_j_new,j],
                                b_out=b_mat[ind_j_new,i],
                                se_exp=se_mat[ind_j_new,j],
                                se_out=se_mat[ind_j_new,i],
                                n = min(n_vec[i],n_vec[j]),
                                rho = rho_ij,
                                random_start = random_start,random_seed=seed)

      obs_graph[i,j] = ItoJ_cML_O_res$BIC_theta
      obs_graph[j,i] = JtoI_cML_O_res$BIC_theta
      obs_graph_se[i,j] = ItoJ_cML_O_res$BIC_se
      obs_graph_se[j,i] = JtoI_cML_O_res$BIC_se
      obs_graph_pval[i,j] = ItoJ_cML_O_res$BIC_p
      obs_graph_pval[j,i] = JtoI_cML_O_res$BIC_p
      k = k + 1

    }
  }

  out = list()
  out$obs_graph = obs_graph
  out$obs_graph_se = obs_graph_se
  out$obs_graph_pval = obs_graph_pval

  return(out)

}

#' @keywords internal
#' @noRd
Graph_Perturb <- function(n_vec,rho_mat,screen_res,
                          num_pert=100,random_start=10,seed=0,trait_vec=NULL){
  set.seed(seed)
  t = 0
  if(is.null(trait_vec)){trait_vec=colnames(screen_res$b_mat)}

  out_list = pbmcapply::pbmclapply(1:num_pert,function(i){
    b_mat_dp = Generate_Perturb(b_mat=screen_res$b_mat,
                                se_mat=screen_res$se_mat,
                                rho_mat=rho_mat,DP_mat_list=screen_res$DP_mat_list);
    Graph_Estimate(b_mat=b_mat_dp,
                   se_mat=screen_res$se_mat,
                   n_vec=n_vec,rho_mat=rho_mat,
                   IJ_snp_list=screen_res$IJ_snp_list,
                   random_start=random_start)},ignore.interactive=TRUE)
  obs_graph_list = lapply(out_list,function(x){x$obs_graph})
  obs_graph_se_list = lapply(out_list,function(x){x$obs_graph_se})
  obs_graph_pval_list = lapply(out_list,function(x){x$obs_graph_pval})

  out = list()
  out$obs_graph_list = obs_graph_list
  out$obs_graph_se_list = obs_graph_se_list
  out$obs_graph_pval_list = obs_graph_pval_list
  out$trait_vec = trait_vec
  return(out)
}

#' @keywords internal
#' @noRd
subset_Graph <- function(dp_list,keep_trait,check=TRUE){
  if(keep_trait=='all'){
    keep_trait_id=1:nrow(dp_list$obs_graph_list[[1]])
  }else{
    keep_trait_id = is.element(dp_list$trait_vec,keep_trait)
  }

  MR_graph_list = dp_list$obs_graph_list
  new_dir_graph_list = lapply(MR_graph_list, function(G_obs) {
    G_dir = Alg1(G_obs, thres = 1e-8)
    G_dir;
  })
  max_eigen_list = lapply(new_dir_graph_list, function(G_dir) {max(abs(eigen(G_dir)$values))})
  rm_dp_id = which(unlist(max_eigen_list)>1)
  if(check & length(rm_dp_id)>0){
    new_dir_graph_list = new_dir_graph_list[-rm_dp_id]
  }
  # new_obs_graph_vec_list = lapply(new_obs_graph_list,function(x){matrix(x,ncol=1)})
  # new_obs_graph_vec_mat = matrix(unlist(new_obs_graph_vec_list),ncol=length(new_obs_graph_vec_list))
  # new_obs_graph_COR = cov(t(new_obs_graph_vec_mat))

  obs_graph_mean = apply(simplify2array(MR_graph_list), 1:2, mean)
  obs_graph_sd = apply(simplify2array(MR_graph_list), 1:2, sd)
  obs_graph_pval = pnorm(-abs(obs_graph_mean/obs_graph_sd))*2

  pval_3dmat = simplify2array(dp_list$obs_graph_pval_list)
  pval_3dmat = pval_3dmat[keep_trait_id,keep_trait_id,]
  dims = dim(pval_3dmat)
  twoDimMat <- matrix(pval_3dmat,prod(dims[1:2]), dims[3])
  twoDimMat[!apply(twoDimMat,1,sum)<1e-100,] -> twoDimMat
  lambda = eigen(cor(t(twoDimMat)))$value
  Me = ceiling(length(lambda) - sum((lambda>1)*(lambda-1)))

  dir_graph_mean = apply(simplify2array(new_dir_graph_list), 1:2, mean)
  dir_graph_sd = apply(simplify2array(new_dir_graph_list), 1:2, sd)
  dir_graph_pval = pnorm(-abs(dir_graph_mean/dir_graph_sd))*2

  colnames(obs_graph_mean)=colnames(obs_graph_sd)=colnames(obs_graph_pval)=dp_list$trait_vec[keep_trait_id]
  colnames(dir_graph_mean)=colnames(dir_graph_sd)=colnames(dir_graph_pval)=dp_list$trait_vec[keep_trait_id]
  rownames(obs_graph_mean)=rownames(obs_graph_sd)=rownames(obs_graph_pval)=dp_list$trait_vec[keep_trait_id]
  rownames(dir_graph_mean)=rownames(dir_graph_sd)=rownames(dir_graph_pval)=dp_list$trait_vec[keep_trait_id]

  eigen_dir = eigen(dir_graph_mean)
  if(any(abs(eigen_dir$values)>1)){warning("Some eigenvalues' absolute values are greater than 1!")}
  dp_list = list()
  dp_list$obs_graph_list = MR_graph_list
  dp_list$dir_graph_list = new_dir_graph_list
  dp_res = list()
  dp_res$MR_graph_Est = obs_graph_mean
  dp_res$MR_graph_sd = obs_graph_sd
  dp_res$MR_graph_pval = obs_graph_pval
  dp_res$Direct_graph_Est = dir_graph_mean
  dp_res$Direct_graph_sd = dir_graph_sd
  dp_res$Direct_graph_pval = dir_graph_pval
  dp_res$Me = Me
  return(list(dp_res=dp_res,dp_list=dp_list))
}



