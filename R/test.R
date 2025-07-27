# p_mat = pnorm(-abs(b_mat/se_mat))*2
# IV_list = list()
# k = 1
# #### select IV #####
# for(i in 1:5){
#   exp_iv = which(p_mat[,i]<5e-8)
#   for(j in (i+1):6){
#     out_iv = which(p_mat[,j]<5e-8)
#     union_iv = union(exp_iv,out_iv)
#
#     data.clump =
#       data.frame(rsid = rownames(p_mat)[union_iv])
#
#     clump.result =
#       ieugwasr::ld_clump(dat = data.clump)
#     IV_list[[k]] = as.character(clump.result$rsid)
#     k = k + 1
#   }
# }

# load('./data/example.rda')
# SP_res = ScreenPair(b_mat = b_mat,se_mat = se_mat,
#                            pval_mat = NULL,n_mat = NULL,
#                            n_vec = n_vec,IV_list = IV_list,R_list = R_list,rho_mat = rho_mat)
# SM_res = ScreenMax(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, pre_screen=SP_res$IJ_snp_list)
# SM_res$DP_mat_list = SP_res$DP_mat_list
# dp_list = Graph_Perturb(screen_res = SM_res,
#                             n_vec = n_vec,
#                             rho_mat = rho_mat,
#                             random_start = 0,
#                             num_pert = 100, seed=1)
#
#
# res = subset_Graph(dp_list=dp_list, keep_trait='all')
# res$dp_res$dir_graph_pval
# B = res$dp_res$dir_graph_pval < 0.05/30
# gg = igraph::graph_from_adjacency_matrix(B, mode = 'directed', weighted=FALSE)
# SA_res = ScreenAug(b_mat = b_mat, se_mat = se_mat, n_vec = n_vec, pre_screen_pair=SP_res$IJ_snp_list,
#                    pre_screen_max=SM_res$IJ_snp_list, graph = gg)
#
# a = MR2G(b_mat = b_mat,se_mat = se_mat,n_vec = n_vec,IV_list = IV_list,R_list = R_list,rho_mat = rho_mat,num_pert_M = 10,num_pert=10,IV_screen = 'ScreenPair')
