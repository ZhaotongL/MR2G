#' @keywords internal
#' @noRd
fn <- function(x,...){
  -loglik(b_t=x,...)
}

#' @keywords internal
#' @noRd
loglik <- function(b_exp,b_out,se_exp,se_out,
                   b_t,theta_t,r_vec_t,rho,t=-qnorm(2.5e-08)){
  l = -1/(2*(1-rho^2)) * sum((
    (b_exp-b_t)^2/se_exp^2 -
      2*rho/(se_exp*se_out)*(b_exp-b_t)*(b_out-theta_t*b_t-r_vec_t) +
      (b_out-theta_t*b_t-r_vec_t)^2/se_out^2
  ))
  logp = sum(log(1 - pnorm(t-b_t/se_exp) + pnorm(-t-b_t/se_exp)))
  return(l-logp)
}

#' @keywords internal
#' @noRd
cML_estimate_O <- function(b_exp,b_out,
                           se_exp,se_out,
                           K,initial_theta = 0,
                           initial_mu = rep(0,length(b_exp)),
                           maxit = 100,rho=0.)
{
  p = length(b_exp)
  ### initialize
  theta = initial_theta
  theta_old = theta-1
  mu_vec = initial_mu

  ###
  ite_ind = 0
  while( (abs(theta_old - theta) > 1e-7) & (ite_ind<maxit))
  {
    theta_old = theta
    ite_ind = ite_ind + 1

    ### first, update v_bg
    if(K>0)
    {
      v_vec = (b_out - theta*mu_vec) - rho*(b_exp-mu_vec)*se_out/se_exp
      v_0 = rep(0,p)
      v_importance = ((b_out - theta*mu_vec)^2 / se_out^2 -
                        2*rho*(b_exp-mu_vec)*(b_out-theta*mu_vec)/(se_exp*se_out)) -
        ((b_out - theta*mu_vec - v_vec)^2 / se_out^2 -
           2*rho*(b_exp-mu_vec)*(b_out-theta*mu_vec-v_vec)/(se_exp*se_out))
      nonzero_bg_ind = sort((order(v_importance,decreasing = T))[1:K])
      v_bg = rep(0,p)
      v_bg[nonzero_bg_ind] = v_vec[nonzero_bg_ind]
    } else{
      v_bg = rep(0,p)
    }

    ### second, update mu_vec
      mu_vec =
        (b_exp/se_exp^2 - rho*(b_out+theta*b_exp-v_bg)/(se_exp*se_out) + theta*(b_out-v_bg)/se_out^2) /
        (1/se_exp^2 - 2*rho*theta/(se_exp*se_out) + theta^2/se_out^2)

    ### third, update theta
    theta =
      sum((b_out - v_bg)*mu_vec / se_out^2 - rho*mu_vec*(b_exp-mu_vec)/(se_exp*se_out)) /
      sum(mu_vec^2 / se_out^2)

  }


  ### one more step for v_bg and mu
  if(K>0)
  {
    nonzero_ind = which(v_bg!=0)
    mu_vec[nonzero_ind] = b_exp[nonzero_ind]
    v_bg[nonzero_ind] = ((b_out - theta*mu_vec) - rho*(b_exp-mu_vec)*se_out/se_exp)[nonzero_ind]
  }

  return(list(theta = theta,
              b_vec = mu_vec,
              r_vec = v_bg))

}

#' @keywords internal
#' @noRd
cML_SdTheta_O <- function(b_exp,b_out,
                          se_exp,se_out,
                          theta,b_vec,r_vec,rho)
{
  t = 0
  nonzero_ind = which(r_vec!=0)
  zero_ind = which(r_vec==0)
  a = sum((b_vec^2/se_out^2)[zero_ind])/(1-rho^2)
  b = 1/(1-rho^2) * ((rho*b_exp-2*rho*b_vec)/(se_exp*se_out) - (b_out-2*theta*b_vec)/se_out^2)
  c = 1/(1-rho^2) * (1/se_exp^2 - 2*rho*theta/(se_exp*se_out) + theta^2/se_out^2)
  f = 1/se_exp * (dnorm(t-b_vec/se_exp) - dnorm(-t-b_vec/se_exp)) ##\partial{p}\partial{bxi}
  p = 1 - pnorm(t-b_vec/se_exp) + pnorm(-t-b_vec/se_exp)
  df = 1/se_exp^2 * ((t+b_vec/se_exp)*dnorm(-t-b_vec/se_exp)+
                       (t-b_vec/se_exp)*dnorm(t-b_vec/se_exp)
  )
  c = c + (-f^2/p^2 + df/p)
  VarTheta = 1/(a - sum((b^2/c)[zero_ind]))


  if(VarTheta<=0)
  {
    return(NaN)
  } else {
    return(sqrt(VarTheta))
  }

}

#' @keywords internal
#' @noRd
cML_estimate_random_O <- function(b_exp, b_out,
                                  se_exp, se_out,
                                  K,random_start = 0,
                                  maxit = 100,rho=0)
{
  p = length(b_exp)
  min_theta_range = min(b_out/b_exp)
  max_theta_range = max(b_out/b_exp)

  theta_v_RandomCandidate = NULL
  sd_v_RandomCandidate = NULL
  l_v_RandomCandidate = NULL
  invalid_RandomCandidate = NULL

  for(random_ind in 1:(1+random_start))
  {
    #      ptm <- proc.time()
    if(random_ind == 1)
    {
      initial_theta = 0
      initial_mu = rep(0,length(b_exp))
      #initial_mu = b_exp
    } else {
      initial_theta = runif(1,min = min_theta_range,max = max_theta_range)
      initial_mu = rnorm(p,mean = b_exp,sd = se_exp)
      #      initial_mu = rep(0,length(b_exp))
    }
    #  print(proc.time() - ptm)
    #    ptm <- proc.time()
    MLE_result =
      cML_estimate_O(b_exp,b_out,
                     se_exp,se_out,
                     K = K,initial_theta = initial_theta,
                     initial_mu = initial_mu,
                     maxit = maxit,rho=rho)

    Neg_l = -loglik(b_exp,b_out,se_exp,se_out,MLE_result$b_vec,MLE_result$theta,MLE_result$r_vec,rho)

    sd_theta = cML_SdTheta_O(b_exp,b_out,
                             se_exp,se_out,
                             MLE_result$theta,
                             MLE_result$b_vec,
                             MLE_result$r_vec,
                             rho = rho)

    theta_v_RandomCandidate = c(theta_v_RandomCandidate,MLE_result$theta)
    sd_v_RandomCandidate = c(sd_v_RandomCandidate,sd_theta)
    l_v_RandomCandidate = c(l_v_RandomCandidate,Neg_l)
    invalid_RandomCandidate = rbind(invalid_RandomCandidate,
                                    as.numeric(MLE_result$r_vec))

    #  print(proc.time() - ptm)
  }
  min_neg_l = which.min(l_v_RandomCandidate)

  theta_est = theta_v_RandomCandidate[min_neg_l]
  sd_est = sd_v_RandomCandidate[min_neg_l]
  l_est = l_v_RandomCandidate[min_neg_l]
  r_est = invalid_RandomCandidate[min_neg_l,]

  return(list(theta = theta_est,
              se = sd_est,
              l = l_est,
              r_est = r_est
  )
  )
}

#' @keywords internal
#' @noRd
mr_cML_O <- function(b_exp,b_out,
                     se_exp,se_out,
                     K_vec = 0:(length(b_exp) - 2),
                     random_start = 0,
                     maxit = 100,
                     random_seed = 0,
                     n, rho=0)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }

  rand_theta = NULL
  rand_sd = NULL
  rand_l = NULL
  invalid_mat = NULL
  for(K_value in K_vec)
  {
    rand_res = cML_estimate_random_O(b_exp = b_exp,
                                     b_out = b_out,
                                     se_exp = se_exp,
                                     se_out = se_out,
                                     K = K_value,
                                     random_start = random_start,
                                     maxit = maxit,
                                     rho = rho)
    rand_theta = c(rand_theta,rand_res$theta)
    rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
  }

  ### get result
  theta_v = rand_theta
  sd_v = rand_sd
  l_v = rand_l

  # cML-MA-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)

  # cML-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_theta = theta_v[min_ind]
  BIC_se = sd_v[min_ind]
  BIC_p = pnorm(-abs(BIC_theta/BIC_se))*2
  BIC_invalid = which(invalid_mat[min_ind,]!=0)


  return(list(
              BIC_theta = BIC_theta,
              BIC_se = BIC_se,
              BIC_p = BIC_p,
              BIC_invalid = BIC_invalid,
              l_vec = l_v,
              BIC_vec = log(n) * K_vec + 2 * l_v)
  )
}

