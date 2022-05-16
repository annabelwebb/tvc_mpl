# control function ---- 
tvc_mpl_control = function(lambda = NULL, kappa = 1/0.6, iter = c(1,1000), n_knots = 8, range = c(0.05,0.95),
                           reg_conv = 1e-6, par_initial = c(0,0,1), line_search = c(1, 1, 1), min_knot = 0){
  
  #check omega's will be between 0 and 1
  if(1/kappa < 0 | 1/kappa > 1){
    stop("kappa value is mis-specified")
  }else{
    kappa = kappa
  }
  
  #set up whether lambda is optimised
  if(is.null(lambda)){
    iter[1] = 1
  }else{
    lambda = lambda
  }
  
  #check initial parameter values
  if(par_initial[3] <= 0){
    par_initial[3] = 0.1
  }else{
    par_initial[3] = par_initial[3]
  }
  
  out = list(lambda = lambda, kappa = kappa, iter = iter, n_knots = n_knots, range = range, 
             reg_conv = reg_conv, par_initial = par_initial, line_search = line_search, min_knot = min_knot)
  return(out)
  
}


# data and parameter set up ------
censor_matrix = function(t1,delta){
  out = matrix(FALSE,nrow=length(t1),ncol=4)
  for(i in 1:length(t1)){
    out[i,(delta[i]+1)]=TRUE
  }
  return(out)
}


model_matrix = function(dat, dat.baseline, dat.left){
  X = as.matrix(dat.baseline[,c(2,8,9,13)])
  Z = as.matrix(dat[,c(14,15)]) #CHANGE THIS LATER
  X_long = as.matrix(dat[,c(2,8,9,13)])
  Z_last = Z[which(dat$last_record ==1),]
  rep = dat.baseline$rep
  last_rec = dat$last_record
  left_pos = dat$left_position
  left_z = as.matrix(dat.left[,c(14,15)]) #CHANGE THIS LATER
  left_x = as.matrix(dat.left[,c(2,8,9,13)])
  
  out = list(X = X, Z = Z, X_long = X_long, Z_last = Z_last,rep = rep, last_rec = last_rec, 
             left_pos = left_pos, left_z = left_z, left_x = left_x)
  return(out)
}


par_init = function(model_mat, n_k = ctrl$n_knots, initial = ctrl$par_initial){
  p = ncol(model_mat$X)
  q = ncol(model_mat$Z)
  m = n_k + 3
  df = nrow(model_mat$X)
  dimen = list(p = p, q = q, m = m, df = df)
  
  beta = matrix(initial[1], nrow = p, ncol = 1)
  gamma = matrix(initial[2], nrow = q, ncol = 1)
  theta = matrix(initial[3], nrow = m, ncol = 1)
  
  out = list(dimen = dimen, beta = beta, gamma = gamma, theta = theta)
  return(out)
  
}

# knots, basis functions and penalty -------

#knots
knots_f = function(times, n_k = ctrl$n_knots, range = ctrl$range, min = 0){
  int_knots = quantile(times,seq(range[1],range[2],length.out=n_k))
  bound_knots = c(min, max(times) + 1e-6)
  return(list(int_knots = int_knots, bound_knots=bound_knots))
}

#M-splines
psi_f = function(events, knts = kn$int_knots, Boundary.knots = kn$bound_knots){
  tmp=mSpline(events, degree = 3, knots = knts, Boundary.knots = Boundary.knots, intercept = F)
  tmp = ifelse(tmp<0, 0, tmp)
  return(tmp)
  
  
}

#I-splines
Psi_f = function(events, knts = kn$int_knots, Boundary.knots = kn$bound_knots){
  tmp=iSpline(events, degree = 3, knots = knts, Boundary.knots = Boundary.knots, intercept = F)
  tmp = ifelse(tmp<0, 0, tmp)
  return(tmp)
}

#sum of I-splines for left and right
NPsi_Diff_f = function(i_long, ev, kn){
  NPsi =  Psi_f(events = ev, knts = kn$int_knots, Boundary.knots = kn$bound_knots)
  
  NPsiDiff = data.frame(cbind(i_long,NPsi)) %>%
    group_by(i_long) %>% mutate_at(vars(-group_cols()),list(~ifelse(row_number()==1,.,. - lag(.)))) %>%
    ungroup(i_long) %>% dplyr::select(-i_long) %>% as.matrix
  
  NPsiDiff[NPsiDiff<0] = 0
  
  return(NPsiDiff)
  
}

#R and Rstar
penalty_f = function(order, IntKnt, bryKnt, dimension, norm = 2){
  if(norm == 2){
    
    ordSp = order
    dgrSp = ordSp - 1
    numIntKnt = length(IntKnt)
    numSp = numIntKnt+ordSp
    
    minTime = min(bryKnt)
    maxTime = max(bryKnt)
    
    R=matrix(0, nrow=numSp, ncol=numSp) 
    xknots = c(rep(minTime, ordSp), IntKnt, rep(maxTime, ordSp))
    for (ii in 1:numSp){
      for (jj in ii:numSp){
        if (jj - ii<ordSp){
          kntset = xknots[xknots>=xknots[jj] & xknots<=xknots[ii+ordSp]]
          kntsum = 0
          for (kk in 1:(length(kntset)-1)){
            kntsum = kntsum + mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, Boundary.knots=bryKnt, 
                                      derivs=dgrSp)[ii]*mSpline(kntset[kk], knots=IntKnt, degree=dgrSp, intercept=T, 
                                                                Boundary.knots=bryKnt,derivs=dgrSp)[jj]*(kntset[kk+1]-kntset[kk])
          }
          R[ii, jj] = kntsum
        }
      }
    }
    R[lower.tri(R, diag = FALSE)] = t(R)[lower.tri(R, diag = FALSE)]
    
    p = dimension$p
    q = dimension$q
    m = numSp
    
    Rstar = rbind(matrix(0,p+q,p+q+m),cbind(matrix(0,m,p+q),R))
    
  }else if(norm == 1){
    m = dimension$m
    
    R = matrix(0, m, m)
    diag(R) = 1
    
    p = dimension$p
    q = dimension$q
    Rstar = rbind(matrix(0,p+q,p+q+m),cbind(matrix(0,m,p+q),R))
    
    
  }
  
  
  
  out = list(R = R, Rstar = Rstar)
  return(out)
  
  
}

# likelihood element updates -------

#compute baseline hazard
h0_f = function(psi, theta){
  h0 = psi%*%theta
  return(h0)
}

#compute Psi star (left and right)

PsiStar_f = function(gamma, i_long, Z, NPsiDiff, tal){
  zGamma = Z%*%gamma
  ezGamma = exp(zGamma)
  m_mat=matrix(rep(1,ncol(NPsiDiff)), nrow=1)
  eZGammaRep = ezGamma %*% m_mat
  
  Psi_star = data.frame(cbind(i_long, eZGammaRep*NPsiDiff)) %>% 
    group_by(i_long) %>% summarise(across(everything(),~sum(.x))) %>% dplyr::select(-i_long) %>% as.matrix
  
  Psi_star_long = NULL
  for(i in 1:length(tal)){
    Psi_star_long = rbind(Psi_star_long, matrix(rep(Psi_star[i,], tal[i]), nrow = tal[i], byrow = TRUE))
  }
  
  out = list(Psi_star = Psi_star, Psi_star_long = Psi_star_long)
  return(out)
}


#compute cumulative hazard function(s)
H_f = function(Psi_Star, beta, theta, X, tal){
  xBeta = X%*%beta
  H0_star = (Psi_Star$Psi_star)%*%theta
  H = H0_star * exp(xBeta)
  
  #long
  H_long = NULL
  for(i in 1:length(tal)){
    H_long = c(H_long, rep(H[i], tal[i]))
  }
  
  out = list(H = H, H_long = H_long)
  return(out)
  
  
}

#compute survival function(s)
S_f = function(H){
  S = exp(-H$H)
  S_long = exp(-H$H_long)
  
  out = list(S = S, S_long = S_long)
  return(out)
}



#compute penalty term info
pen_term_f = function(penalty, theta, lambda){
  R = penalty$R
  Rtheta = R%*%theta
  thetaRtheta = t(theta)%*%Rtheta
  
  #print(c(lambda, Rtheta))
  
  TwoLRtheta = lambda*2*Rtheta
  
  penalise_ll = lambda*thetaRtheta
  
  out = list(Rtheta = Rtheta, thetaRtheta = thetaRtheta, TwoLRtheta = TwoLRtheta, penalise_ll = penalise_ll)
  return(out)
  
}

pen_lik_f = function(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, pen_term, model_mat, beta, gamma, theta, censor){
  X = model_mat$X
  Z = model_mat$Z
  
  #event
  xBeta = X%*%beta
  z_last = model_mat$Z_last
  zGamma_last = z_last%*%gamma
  event_ll = sum((log(h0)[censor[,2]] + xBeta[censor[,2]] + zGamma_last[censor[,2]] + log(S_r$S[censor[,2]])))
  
  #right
  right_ll = sum(log(S_r$S)[censor[,1]])
  
  #left
  S_l_diff = (1 - S_r$S)[censor[,3]]
  left_ll = sum(log(S_l_diff))
  
  #interval
  S_diff = (S_l$S[censor[,4]] - S_r$S[censor[,4]])
  S_diff[S_diff<=0] = 1e-4
  interval_ll = sum(log(S_diff))
  
  #print(c(event_ll, right_ll))
  
  #compute
  pen_lik = event_ll + right_ll + left_ll + interval_ll - pen_term$penalise_ll
  return(pen_lik)
}

# score and hessian computation -------

beta_score = function(X, S_r, S_l, H_r, H_l, censor){
  beta_score = t(X[censor[,2],]) %*% (1 - H_r$H[censor[,2]]) - 
    t(X[censor[,1],]) %*% (H_r$H[censor[,1]]) +
    t(matrix(X[censor[,3],], ncol = ncol(X))) %*% as.matrix((S_r$S[censor[,3]] * H_r$H[censor[,3],] / (1 - S_r$S[censor[,3]]))) - 
    t(X[censor[,4],]) %*% ((S_l$S[censor[,4]] * H_l$H[censor[,4]] - S_r$S[censor[,4]] * H_r$H[censor[,4]]) / (S_l$S[censor[,4]] - S_r$S[censor[,4]]) )
  
  
  return(beta_score)
}

beta_hessian = function(X, S_r, S_l, H_r, H_l, censor){
  beta_neghess = t(X[censor[,2],]) %*% diag(as.numeric(H_r$H[censor[,2]])) %*% (X[censor[,2],]) +
    t(X[censor[,1],]) %*% diag(as.numeric(H_r$H[censor[,1]])) %*% (X[censor[,1],]) -
    t(matrix(X[censor[,3],], ncol = ncol(X))) %*% diag(as.numeric((S_r$S[censor[,3]] * H_r$H[censor[,3]] / (1 - S_r$S[censor[,3]]))), nrow = sum(censor[,3])) %*% (X[censor[,3],]) +
    t(matrix(X[censor[,3],], ncol = ncol(X))) %*% diag(as.numeric((S_r$S[censor[,3]] * (H_r$H^2)[censor[,3]] / ((1 - S_r$S[censor[,3]])^2))), nrow = sum(censor[,3])) %*% (X[censor[,3],]) +
    t(X[censor[,4],]) %*% diag(as.numeric((S_l$S[censor[,4]] * H_l$H[censor[,4]] - S_r$S[censor[,4]] * H_r$H[censor[,4]]) / (S_l$S[censor[,4]] - S_r$S[censor[,4]]) )) %*% (X[censor[,4],]) +
    t(X[censor[,4],]) %*% diag(as.numeric((S_l$S[censor[,4]] * S_r$S[censor[,4]] * (H_l$H[censor[,4]] - H_r$H[censor[,4]])^2) / (S_l$S[censor[,4]] - S_r$S[censor[,4]])^2)) %*% (X[censor[,4],])
  
  
  return(beta_neghess)
}

A_star_f = function(X_long, Z, beta, gamma, theta,NPsiDiff){
  eXBeta = exp(X_long%*%beta)
  zGamma = Z%*%gamma
  zGamma1 = gamma[1] %*% Z[,1]
  zGamma2 = gamma[2] %*% Z[,2]
  
  ezGamma = as.numeric(exp(zGamma))
  ezGamma1 = as.numeric(exp(zGamma1))
  ezGamma2 = as.numeric(exp(zGamma2))
  
  H0_diff = NPsiDiff %*% theta #Nx1
  A_star = H0_diff * eXBeta * ezGamma
  A_star1 = H0_diff * eXBeta * ezGamma1
  A_star2 = H0_diff * eXBeta * ezGamma2
  out = list(A_star = A_star, A_star1 = A_star1, A_star2 = A_star2)
  
  return(out)
  
}

gamma_score = function(model_mat, censor_long, A_r, A_l, S_r, S_l, last){
  Z = model_mat$Z
  
  oneMinusS = (1 - S_r$S_long)+ 1e-3
  SDiff = (S_l$S_long - S_r$S_long) + 1e-3
  
  d_1 = as.numeric(censor_long[,1] + censor_long[,2] - censor_long[,3]*S_r$S_long/(oneMinusS) - 
                     censor_long[,4]*S_r$S_long/(SDiff))
  d_2 =  censor_long[,4]*as.numeric(S_l$S_long/(SDiff))
  
  A_mat_r = cbind(A_r$A_star1, A_r$A_star2)
  A_mat_l = cbind(A_l$A_star1, A_l$A_star2)
  
  gamma_score = t(Z) %*% (last*as.numeric(censor_long[,2])) - 
    (t(Z) %*% (d_1*A_mat_r))[,1] -
    (t(Z) %*% (d_2*A_mat_l))[,1]
  
  return(gamma_score)
  
}


gamma_hessian = function(model_mat, censor_long, A_r, A_l, S_r, S_l){
  Z = model_mat$Z
  
  oneMinusS = (1 - S_r$S_long) + 1e-3
  SDiff = (S_l$S_long-S_r$S_long) + 1e-3
  
  d_1 = censor_long[,1] + censor_long[,2] - censor_long[,3]*S_r$S_long/(oneMinusS) - censor_long[,4]*S_r$S_long/(SDiff)
  d_2 = censor_long[,4]*S_l$S_long/(SDiff)
  d_3 = censor_long[,3]*S_r$S_long/((oneMinusS)^2)
  d_4 = censor_long[,4]*S_l$S_long*S_r$S_long/((SDiff)^2)
  
  
  gamma_neghess = t(Z) %*% diag(as.numeric(d_1*A_r$A_star)) %*% Z +
    t(Z) %*% diag(as.numeric(d_2*A_l$A_star)) %*% Z +
    t(Z) %*% diag(as.numeric(d_3* A_r$A_star1 * A_r$A_star2)) %*% Z +
    t(Z) %*% diag(as.numeric(d_4 * (A_l$A_star1 - A_r$A_star1) * (A_l$A_star2 - A_r$A_star2))) %*% Z
  
  return(gamma_neghess)
  
}

theta_score = function(psi, h0, S_r, S_l, PsiStar_r, PsiStar_l, X, beta, p_t, censor){
  eXBeta = exp(X%*%beta)
  
  delta = any(censor[,2])
  
  theta_score_pos = delta*(t(psi[censor[,2],]) %*% (1/h0[censor[,2],])) + 
    t(matrix(PsiStar_r$Psi_star[censor[,3],], ncol = ncol(PsiStar_r$Psi_star))) %*% (S_r$S[censor[,3]] * eXBeta[censor[,3]] / (1 - S_r$S[censor[,3]])) +
    t(PsiStar_r$Psi_star[censor[,4],]) %*% (S_r$S[censor[,4]] * eXBeta[censor[,4]] / (S_l$S[censor[,4]] - S_r$S[censor[,4]]))
  
  theta_score_neg = t(PsiStar_r$Psi_star[censor[,2],]) %*% (eXBeta[censor[,2]]) +
    t(PsiStar_r$Psi_star[censor[,1],]) %*% (eXBeta[censor[,1]]) +
    t(PsiStar_l$Psi_star[censor[,4],]) %*% (S_l$S[censor[,4]] * eXBeta[censor[,4]] / (S_l$S[censor[,4]] - S_r$S[censor[,4]]))
  
  
  TwoLRtheta = p_t$TwoLRtheta
  theta_score_neg = theta_score_neg + TwoLRtheta*(TwoLRtheta>0) + 1e-7
  theta_score_pos = theta_score_pos - TwoLRtheta*(TwoLRtheta<0) + 1e-7
  
  theta_score = theta_score_pos - theta_score_neg
  out = list(theta_score = theta_score, theta_score_pos = theta_score_pos, theta_score_neg = theta_score_neg)
  return(out)
  
}

theta_hessian = function(psi, h0, PsiStar_r, PsiStar_l, X, beta, S_r, S_l, censor){
  
  delta = any(censor[,2])
  eXBeta = exp(X%*%beta)
  d_1 = delta * as.numeric((1/(h0^2)))
  d_2 = as.numeric(censor[,3] * (eXBeta^2) * S_r$S/((1-S_r$S)^2))
  d_3 = as.numeric(censor[,4] * (eXBeta^2) * S_l$S * S_r$S/((S_l$S - S_r$S)^2))
  
  theta_neghess = t(psi) %*% diag(d_1) %*% psi +
    t(PsiStar_r$Psi_star) %*% diag(d_2) %*% PsiStar_r$Psi_star +
    t(PsiStar_l$Psi_star - PsiStar_r$Psi_star) %*% diag(d_3) %*% (PsiStar_l$Psi_star - PsiStar_r$Psi_star)
  
  return(theta_neghess)
  
}


P_star_f = function(beta, gamma, Z, X_long, NPsiDiff){
  eXBeta = exp(X_long%*%beta)
  zGamma = Z%*%gamma
  ezGamma = as.numeric(exp(zGamma))
  
  NDiff_eZGamma = NPsiDiff * ezGamma
  P_star = NDiff_eZGamma * as.numeric(eXBeta)
  return(P_star)
  
}

beta_gamma_hessian = function(censor_long, beta, model_mat, A_r, A_l, S_r, S_l, H_r, H_l){
  Z = model_mat$Z #Nxq
  X = model_mat$X_long #Nxp
  
  d_1 = censor_long[,1] + censor_long[,2] - censor_long[,3]*(S_r$S_long - S_r$S_long*H_r$H_long - S_r$S_long^2)/((1 - S_r$S_long)^2) - 
    censor_long[,4]*(S_r$S_long/(S_l$S_long-S_r$S_long) + S_l$S_long*S_r$S_long*(H_l$H_long-H_r$H_long)/((S_l$S_long-S_r$S_long)^2) )
  
  d_2 = censor_long[,4] * (S_l$S_long/(S_l$S_long-S_r$S_long) + S_l$S_long*S_r$S_long*(H_l$H_long-H_r$H_long)/((S_l$S_long-S_r$S_long)^2) )
  
  A_mat_r = cbind(A_r$A_star1, A_r$A_star2)
  
  beta_gamma = t(Z) %*% diag(as.numeric((A_r$A_star * d_1))) %*% X +
    t(Z) %*% diag(as.numeric((A_l$A_star * d_2))) %*% X
  
  return(beta_gamma)
  
}

beta_theta_hessian = function(X, S_r, S_l, H_r, H_l, PsiStar_r, PsiStar_l, beta, censor){
  eXBeta = exp(X%*%beta)
  d_1 = as.numeric(eXBeta*(censor[,1] + censor[,2] - censor[,3]*S_r$S/(1-S_r$S) + 
                             censor[,3]*S_r$S*H_r$H/((1-S_r$S)^2) - censor[,4]*S_r$S/(S_l$S - S_r$S) -
                             censor[,4] * S_l$S*S_r$S*(H_l$H - H_r$H)/((S_l$S - S_r$S)^2)))
  
  
  d_2 = as.numeric(eXBeta*(censor[,4]*S_l$S/(S_l$S - S_r$S) + 
                             censor[,4] * S_l$S*S_r$S*(H_l$H - H_r$H)/((S_l$S - S_r$S)^2)))
  
  
  beta_theta = t(X) %*% (d_1 * PsiStar_r$Psi_star) +
    t(X) %*% (d_2 * PsiStar_l$Psi_star)
  
  return(beta_theta)
  
}

gamma_theta_hessian = function(P_r, P_l, A_r, A_l, censor_long, S_r, S_l, PsiStar_r, PsiStar_l, X_long, Z, beta){
  
  eXBeta = exp(X_long%*%beta)
  
  d_1 = as.numeric(censor_long[,2] + censor_long[,1] - censor_long[,3]*S_r$S_long/(1-S_r$S_long) - 
                     censor_long[,4]*S_r$S_long/(S_l$S_long-S_r$S_long))
  
  d_2 = as.numeric(censor_long[,4]*S_l$S_long/(S_l$S_long-S_r$S_long))
  
  d_3 = as.numeric(censor_long[,3]*eXBeta*S_r$S_long*A_r$A_star/((1-S_r$S_long)^2))
  
  d_4 = as.numeric(censor_long[,4]*eXBeta*S_l$S_long*S_r$S_long*(A_l$A_star - A_r$A_star)/((S_l$S_long-S_r$S_long)^2))
  
  gamma_theta = t(Z) %*% (P_r * d_1) +
    t(Z) %*% (P_l * d_2) + 
    t(Z) %*% (PsiStar_r$Psi_star_long * d_3) + 
    t(Z) %*% ((PsiStar_l$Psi_star_long - PsiStar_r$Psi_star_long) * d_4)
  
  return(gamma_theta)
  
}

# parameter updates --------

#update beta
update_beta = function(beta_old, censor, model_mat, S_r, S_l, H_r, H_l){
  X = model_mat$X
  beta_score = beta_score(X, S_r, S_l, H_r, H_l, censor)
  beta_neghess = beta_hessian(X, S_r, S_l, H_r, H_l, censor)
  
  step = solve(beta_neghess)%*%beta_score
  
  beta = beta_old + step
  
  out = list(beta = beta, step = step)
  
  return(out)
  
}

update_gamma = function(gamma_old, model_mat, A_r, A_l, S_r, S_l, censor_long, last){
  
  gamma_score = gamma_score(model_mat, censor_long, A_r, A_l, S_r, S_l, last)
  gamma_neghess = gamma_hessian(model_mat, censor_long, A_r, A_l, S_r, S_l)
  #print(gamma_neghess)
  step = solve(gamma_neghess)%*%gamma_score
  
  gamma = gamma_old + step
  
  #print(c(gamma_score, solve(gamma_neghess)))
  
  out = list(gamma = gamma, step = step, gs = gamma_score, gh = gamma_neghess)
  
  return(out)
  
}

update_theta = function(theta_old, X, psi, h0, PsiStar_r, PsiStar_l, beta, p_t, S_r, S_l, censor, ac_con, m){
  
  th_s = theta_score(psi, h0, S_r, S_l, PsiStar_r, PsiStar_l, X, beta, p_t, censor)
  th_neg = th_s$theta_score_neg
  th_pos = th_s$theta_score_pos
  theta_score = th_s$theta_score
  
  D_matrix = matrix(0, nrow = m, ncol = 1)
  D_matrix[ac_con] = theta_old[ac_con]/th_neg[ac_con]
  
  step = D_matrix*(theta_score)
  
  theta = theta_old + step
  #print(step)
  
  out = list(theta = theta, step = step, score = theta_score, neg_score = th_neg, pos_score = th_pos)
  
  return(out)
  
}


zTA_i_f = function(Z, A_r, A_l, n, i_long){
  zTA_r = matrix(0, nrow = n, ncol = ncol(Z))
  zTA_l = matrix(0, nrow = n, ncol = ncol(Z))
  
  for(i in 1:n){
    z_i = Z[which(i_long == i),]
    A_r_i = A_r$A_star[which(i_long == i),]
    zTA_r[i,] = A_r_i %*% (z_i)
    A_l_i = A_l$A_star[which(i_long == i),]
    zTA_l[i,] = A_l_i %*% (z_i)
    
  }
  
  out = list(zTA_r = zTA_r, zTA_l = zTA_l)
  return(out)
  
}





B_matrix = function(Z, gamma, NPsiDiff, rep, theta){
  
  #for each i, get row vector of H0 diffs across each time 
  zGamma = Z%*%gamma
  ezGamma = as.numeric(exp(zGamma))
  
  H0_diff = NPsiDiff %*% theta
  H0_diff = H0_diff * ezGamma
  
  B = matrix(0, ncol = nrow(Z), nrow = length(rep))
  for(r in 1:length(rep)){
    if(r == 1){
      B_i = H0_diff[1:(rep[r]),]
      B[r,1:(rep[r])] = B_i
      pos = rep[r]
    }else{
      B_i = H0_diff[(pos + 1):(pos + rep[r]),]
      B[r,(pos + 1):(pos + rep[r])] = B_i
      pos = pos + rep[r]
    }
  }
  
  return(B)
  
}


#D matrix

D_matrix = function(Z, gamma, NPsiDiff, theta){
  
  zGamma = Z%*%gamma
  ezGamma = as.numeric(exp(zGamma))
  
  H0_diff = NPsiDiff %*% theta
  D = H0_diff * ezGamma
  
  return(D)
  
}


#score
gamma_score_jun = function(B_r, B_l, censor, S_r, S_l, model_mat, beta, censor_long){
  
  X = model_mat$X
  Z = model_mat$Z
  eXBeta = exp(X%*%beta)
  epsilon = model_mat$last_rec*censor_long[,2]
  
  d1 = matrix(censor[,2]*eXBeta + censor[,1]*eXBeta - censor[,3]*eXBeta*S_r$S/(1-S_r$S) )
  d2 = matrix(censor[,4]*eXBeta*S_l$S/(S_l$S-S_r$S))
  d3 = matrix(censor[,4]*eXBeta*S_r$S/(S_l$S - S_r$S))
  
  mat = t(B_r) %*% d1 + t(B_l) %*% d2 - t(B_r) %*% d3
  
  gamma_score = t(Z) %*% (epsilon - mat)
  #print(c(t(Z) %*% epsilon, t(Z) %*% (t(B$B_right) %*% d1 + t(B$B_left) %*% d2)))
  return(gamma_score)
  
}

#hessian

gamma_hess_jun = function(D_r, D_l, B_r, B_l, censor_long, mod_mat, beta, S_r, S_l){
  eXBeta = exp(mod_mat$X_long%*%beta)
  Z = mod_mat$Z
  
  d1 = diag(as.numeric((censor_long[,1] + censor_long[,2]) * eXBeta * D_r))
  d2 = diag(as.numeric(censor_long[,3] * D_r * eXBeta * S_r$S_long/(1 - S_r$S_long)))
  d3 = diag(as.numeric(censor_long[,4] * S_l$S_long * eXBeta * D_r/(S_l$S_long - S_r$S_long)))
  d4 = diag(as.numeric(censor_long[,4] * D_r * eXBeta * S_r$S_long / (S_l$S_long - S_r$S_long)))
  
  b1 = censor_long[,3] * eXBeta * eXBeta * S_r$S_long/(1 - S_r$S_long)
  b2 = censor_long[,4] * eXBeta * eXBeta * S_l$S_long * S_r$S_long / (S_l$S_long - S_r$S_long)
  
  BTB_r = t(B_r) %*% B_r
  BTB_rl = t(B_r- B_l) %*% (B_r - B_l)
  
  for(t in 1:nrow(BTB_r)){
    BTB_r[,t] = BTB_r[,t] * b1[t]
    BTB_rl[,t] = BTB_rl[,t] * b2[t]
    
    
  }
  
  Gamma = d1 + BTB_r + BTB_rl + d3
  
  g_h = t(Z) %*% Gamma %*% Z
  return(g_h)
  
  
  
}


update_gamma_jun = function(gamma_old, model_mat, S_r, S_l, beta, censor_long, B_r, B_l, D_r, D_l, censor){
  
  gamma_score = gamma_score_jun(B_r, B_l, censor, S_r, S_l, model_mat, beta, censor_long)
  gamma_neghess = gamma_hess_jun(D_r, D_l, B_r, B_l, censor_long, model_mat, beta, S_r, S_l)
  
  step = solve(gamma_neghess)%*%gamma_score
  
  gamma = gamma_old + step
  
  #print(c(gamma_score, solve(gamma_neghess)))
  
  out = list(gamma = gamma, step = step, gs = gamma_score, gh = gamma_neghess)
  
  return(out)
  
}

