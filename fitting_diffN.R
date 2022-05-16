tvc_fit = function(dat, dat.baseline, dat.left, ctrl){
  
  mod_mat = model_matrix(dat,dat.baseline, dat.left)
  
  time1 = dat.baseline$TL
  time2 = dat.baseline$TR
  
  censor = censor_matrix(time1, dat.baseline$delta)
  
  time1_long = as.numeric(dat$TL_long)
  time2_long = as.numeric(dat$TR_long)
  
  censor_long = censor_matrix(time1_long, dat$delta_long)
  
  #get knots
  times = c(time1[censor[,1]],time1[censor[,2]],
            time1[censor[,4]],time2[censor[,4]],time2[censor[,3]])
  
  
  kn = knots_f(times = times, n_k = ctrl$n_knots, range = ctrl$range, min = 0)
  init = par_init(mod_mat, initial = ctrl$par_initial)
  beta = init$beta
  gamma = init$gamma
  theta = init$theta
  df = init$dimen$df
  kappa = ctrl$kappa
  
  m = init$dimen$m
  p = init$dimen$p
  q = init$dimen$q
  
  #get basis functions
  times_psi = dat.baseline$TR
  times_psi[censor[,1]] = dat.baseline$TL[censor[,1]]
  psi = psi_f(times_psi, knts = kn$int_knots, Boundary.knots = kn$bound_knots)
  
  
  NPsiDiff_r = NPsi_Diff_f(dat$i_long, ev = dat$end, kn = kn)
  NPsiDiff_l = NPsi_Diff_f(dat.left$i_long, ev = dat.left$end, kn = kn)
  
  
  
  lambda = ctrl$lambda
  penalty = penalty_f(3, kn$int_knots, kn$bound_knots, init$dimen, norm = 2)
  warn = 0
  
  for(it in 1:(ctrl$iter[1])){ #start outer loop
    ac_con = rep(TRUE, m)
    for(iter in 1:(ctrl$iter[2])){ #start inner loop
      
      #initial values
      h0 = h0_f(psi, theta)
      
      PsiStar_r = PsiStar_f(gamma, dat$i_long, mod_mat$Z, NPsiDiff_r, dat.baseline$rep)
      PsiStar_l = PsiStar_f(gamma, dat.left$i_long, mod_mat$left_z, NPsiDiff_l, dat.baseline$rep)
      
      H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
      H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
      
      
      S_r = S_f(H_r)
      S_l = S_f(H_l)
      
      p_t = pen_term_f(penalty, theta, lambda)
      
      log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
      
      #update beta
      loglik_OLD = log_lik
      beta_OLD = beta
      beta_new = update_beta(beta_old = beta_OLD, censor, mod_mat, S_r, S_l, H_r, H_l)
      beta = beta_new$beta
      
      
      #new likelihood
      H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
      H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
      
      
      S_r = S_f(H_r)
      S_l = S_f(H_l)
      
      log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
      
      #step size
      if((log_lik<loglik_OLD) & ctrl$line_search[1] == 1){
        i = 0
        omega = 1/kappa
        while(log_lik < loglik_OLD){
          H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
          H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
          
          
          S_r = S_f(H_r)
          S_l = S_f(H_l)
          
          log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
          
          if(omega>=1e-2){
            omega = omega/kappa
          }else if(omega<1e-2 & omega>=1e-5){
            omega = omega*(5e-2)
          }else if(omega<1e-5){
            omega = omega*(1e-5)
          }
          
          i = i+1
          if(i>500){break}
        }
        
        
      }
      
      
      #update gamma
      B_r = B_matrix(mod_mat$Z, gamma, NPsiDiff_r, dat.baseline$rep, theta)
      B_l = B_matrix(mod_mat$Z, gamma, NPsiDiff_l, dat.baseline$rep, theta)
      
      
      D_r = D_matrix(mod_mat$Z, gamma, NPsiDiff_r, theta)
      D_l = D_matrix(mod_mat$Z, gamma, NPsiDiff_l, theta)
      
      A_s_r = A_star_f(mod_mat$X_long, mod_mat$Z, beta, gamma, theta, NPsiDiff_r)
      A_s_l = A_star_f(mod_mat$left_x, mod_mat$left_z, beta, gamma, theta, NPsiDiff_l)
      
      gamma_OLD = gamma
      
      gamma_new = update_gamma(gamma_OLD, mod_mat, A_s_r, A_s_l, S_r, S_l, censor_long, mod_mat$last_rec)
      
      gamma = gamma_new$gamma
      
      
      #new likelihood
      PsiStar_r = PsiStar_f(gamma, dat$i_long, mod_mat$Z, NPsiDiff_r, dat.baseline$rep)
      PsiStar_l = PsiStar_f(gamma, dat.left$i_long, mod_mat$left_z, NPsiDiff_l, dat.baseline$rep)
      
      H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
      H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
      
      
      S_r = S_f(H_r)
      S_l = S_f(H_l)
      
      log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
      
      
      #step size
      if((log_lik < loglik_OLD) & ctrl$line_search[2] == 1){
        #print(c(loglik_OLD - log_lik))
        i = 0
        omega = 1/kappa
        
        while(log_lik < loglik_OLD){
          #print(omega)
          gamma = gamma_OLD + omega*gamma_new$step
          PsiStar_r = PsiStar_f(gamma, dat$i_long, mod_mat$Z, NPsiDiff_r, dat.baseline$rep)
          PsiStar_l = PsiStar_f(gamma, dat.left$i_long, mod_mat$left_z, NPsiDiff_l, dat.baseline$rep)
          
          H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
          H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
          
          S_r = S_f(H_r)
          S_l = S_f(H_l)
          
          log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
          
          if(omega>=1e-2){
            omega = omega/kappa
          }else if(omega<1e-2 & omega>=1e-5){
            omega = omega*(5e-2)
          }else if(omega<1e-5){
            omega = omega*(1e-5)
          }
          
          
          i = i+1
          if(i>500){break}
        }
      }
      
      
      
      #update theta
      loglik_OLD = log_lik
      
      
      theta_OLD = theta
      
      theta_new = update_theta(theta_OLD, mod_mat$X, psi, h0, PsiStar_r, PsiStar_l, beta, p_t, S_r, S_l, censor,  ac_con, m)
      theta = theta_new$theta
      
      
      #new likelihood
      
      h0 = h0_f(psi, theta)
      H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
      H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
      
      
      S_r = S_f(H_r)
      S_l = S_f(H_l)
      
      p_t = pen_term_f(penalty, theta, lambda)
      
      log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
      
      #step size
      if((log_lik < loglik_OLD) & ctrl$line_search[3] == 1){
        i = 0
        omega = 1/kappa
        while(log_lik < loglik_OLD){
          theta = theta_OLD + omega*theta_new$step
          h0 = h0_f(psi, theta)
          H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
          H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
          
          
          S_r = S_f(H_r)
          S_l = S_f(H_l)
          
          p_t = pen_term_f(penalty, theta, lambda)
          
          log_lik = pen_lik_f(h0, PsiStar_r, PsiStar_l, H_r, H_l, S_r, S_l, p_t, mod_mat, beta, gamma, theta, censor)
          
          if(omega>=1e-2){
            omega = omega/kappa
          }else if(omega<1e-2 & omega>=1e-5){
            omega = omega*(5e-2)
          }else if(omega<1e-5){
            omega = omega*(1e-5)
          }
          
          i = i+1
          if(i>500){break}
        }
        #print(omega)
      }
      
      for(u in 1:m){
        if(theta[u] < 0.01 & theta_new$score[u] < (-0.1)){
          ac_con[u] = FALSE
          theta[u] < 1e-4
          
        }
      }
      
      A_s_r = A_star_f(mod_mat$X_long, mod_mat$Z, beta, gamma, theta, NPsiDiff_r)
      A_s_l = A_star_f(mod_mat$left_x, mod_mat$left_z, beta, gamma, theta, NPsiDiff_l)
      
      beta_gradient = beta_score(mod_mat$X, S_r, S_l, H_r, H_l, censor)
      gamma_gradient = gamma_score(mod_mat, censor_long, A_s_r, A_s_l, S_r, S_l, dat$last_record)
      
      
      parameters = list(beta = beta, gamma = gamma, theta = theta)
      theta_info = list(th_s = theta_new$score, constraints = which(ac_con == FALSE))
      v = seq(0, max(times), length.out = 1000)
      psi_est = psi_f(v, knts = kn$int_knots, Boundary.knots = kn$bound_knots)
      h0_est = psi_est%*%theta
      Psi_est = Psi_f(v, knts = kn$int_knots, Boundary.knots = kn$bound_knots)
      H0_est = Psi_est%*%theta
      S0_est = exp(-H0_est)
      func_est = list(h0_est = h0_est, H0_est = H0_est, S0_est = S0_est, v = v)
      grad = c(beta_gradient, gamma_gradient)
      print(c(iter, beta, gamma))
      
      #print(c(iter, (theta_new$neg_score[ac_con] - theta_new$pos_score[ac_con])))
      
      #print(c(iter, abs(beta-beta_OLD), abs(gamma-gamma_OLD), abs(theta-theta_OLD), which(ac_con == FALSE)))
      
      if((all(abs(c(beta-beta_OLD, gamma-gamma_OLD, 
                    (theta_new$neg_score[ac_con] - theta_new$pos_score[ac_con]))) < 10*ctrl$reg_conv))){
        reg.conv = 1
        break
      }else{
        reg.conv = 0
      }
      
      
      
      
    } #end inner loop
    #hessian matrix
    #print(c(iter, beta, gamma))
    Hess=HRinv=matrix(0,p+q+m,p+q+m)
    
    PsiStar_r = PsiStar_f(gamma, dat$i_long, mod_mat$Z, NPsiDiff_r, dat.baseline$rep)
    PsiStar_l = PsiStar_f(gamma, dat.left$i_long, mod_mat$left_z, NPsiDiff_l, dat.baseline$rep)
    
    H_r = H_f(PsiStar_r, beta, theta, mod_mat$X, dat.baseline$rep)
    H_l = H_f(PsiStar_l, beta, theta, mod_mat$X, dat.baseline$rep)
    
    S_r = S_f(H_r)
    S_l = S_f(H_l)
    
    A_s_r = A_star_f(mod_mat$X_long, mod_mat$Z, beta, gamma, theta, NPsiDiff_r)
    A_s_l = A_star_f(mod_mat$left_x, mod_mat$left_z, beta, gamma, theta, NPsiDiff_l)
    
    h0 = h0_f(psi, theta)
    
    P_s_r = P_star_f(beta, gamma, mod_mat$Z, mod_mat$X_long, NPsiDiff_r)
    P_s_l = P_star_f(beta, gamma, mod_mat$Z, mod_mat$X_long, NPsiDiff_l)
    
    #beta beta
    Hess[1:p,1:p] = beta_hessian(mod_mat$X, S_r, S_l, H_r, H_l, censor)
    
    #beta gamma
    Hess[1:p,(p+1):(p+q)] = beta_gamma_hessian(censor_long, beta, mod_mat, A_s_r, A_s_l, S_r, S_l, H_r, H_l)
    Hess[(p+1):(p+q),1:p] = t(Hess[1:p,(p+1):(p+q)])
    
    #beta theta
    Hess[1:p,(p+q+1):(p+q+m)] = beta_theta_hessian(mod_mat$X, S_r, S_l, H_r, H_l, PsiStar_r, PsiStar_l, beta, censor)
    Hess[(p+q+1):(p+q+m),1:p]=t(Hess[1:p,(p+q+1):(p+q+m)])
    
    #gamma gamma
    Hess[(p+1):(p+q),(p+1):(p+q)] = gamma_hessian(mod_mat, censor_long, A_s_r, A_s_l, S_r, S_l)
    
    #gamma theta
    Hess[(p+1):(p+q),(p+q+1):(p+q+m)] = gamma_theta_hessian(P_s_r, P_s_l, A_s_r, A_s_l, censor_long, S_r, S_l, PsiStar_r, PsiStar_l, mod_mat$X_long, mod_mat$Z, beta)
    Hess[(p+q+1):(p+q+m),(p+1):(p+q)]=t(Hess[(p+1):(p+q),(p+q+1):(p+q+m)])
    
    #theta theta
    Hess[(p+q+1):(p+q+m),(p+q+1):(p+q+m)] = theta_hessian(psi, h0, PsiStar_r, PsiStar_l, mod_mat$X, beta, S_r, S_l, censor)
    
    
    #update lambda
    
    lambda_old   = lambda
    df_old       = df
    sigma2_old   = 1/(2*lambda_old)
    pos=c(rep(TRUE,p+q+m))
    pos[p+q+which(!ac_con)]=FALSE
    
    R = penalty$R
    Q_mat = as.numeric(1/sigma2_old)*penalty$Rstar
    
    
    if(ctrl$iter[1]==1){
      lambda = ctrl$lambda
    }else{
      HRinv = solve(-Hess+Q_mat)
      df           = m-sum(diag(HRinv%*%Q_mat))
      #print(c(iter, df))
      
      if(df<0){
        
        stop("Negative variance in estimation of smoothing parameter. Number of knots may be mis-specified.")
      }
      sigma2       = as.numeric(t(theta)%*%R%*%theta/df )
      
      #as.numeric(t(theta[ac_con])%*%R[ac_con, ac_con]%*%theta[ac_con]/df )
      
      lambda       = 1/(2*sigma2)
      print(lambda)
    }
    
    
    
    #print(abs(df-df_old))
    
    
    
    
    if((abs(df-df_old)<5)){
      sm.conv=1
      #print(c("converged"))
      break
    }else{
      sm.conv=0
    }
    
    
    
    
    
  } #end outer loop
  
  
  
  p_t = pen_term_f(penalty, theta, lambda)
  
  Rstar = penalty$Rstar # does this need 1/sigma2?
  
  
  
  #try making Q matrix
  Q = matrix(NA,nrow(mod_mat$X),(p+q+m))
  
  zTA = zTA_i_f(mod_mat$Z, A_s_r, A_s_l, nrow(dat.baseline), dat$i_long)
  eXBeta = exp(mod_mat$X%*%beta)
  
  Q[censor[,2],(1:p)] = mod_mat$X[censor[,2],] * (1 - H_r$H[censor[,2]])
  Q[censor[,1],(1:p)] = mod_mat$X[censor[,1],] * (- H_r$H[censor[,1]])
  Q[censor[,3],(1:p)] = mod_mat$X[censor[,3],] * (S_r$S[censor[,3]] * H_r$H[censor[,3],] / (1 - S_r$S[censor[,3]]))
  Q[censor[,4],(1:p)] = mod_mat$X[censor[,4],] * (-((S_l$S[censor[,4]] * H_l$H[censor[,4]] - S_r$S[censor[,4]] * H_r$H[censor[,4]]) / (S_l$S[censor[,4]] - S_r$S[censor[,4]]) ))
  
  Q[censor[,2],(p+1):(p+q)] = mod_mat$Z_last[censor[,2]] - zTA$zTA_r[censor[,2]]
  Q[censor[,1],(p+1):(p+q)] = - zTA$zTA_r[censor[,1]]
  Q[censor[,3],(p+1):(p+q)] = zTA$zTA_r[censor[,3]] * (S_r$S/(1-S_r$S))[censor[,3]]
  Q[censor[,4],(p+1):(p+q)] = - ((as.numeric(S_l$S) * zTA$zTA_l - as.numeric(S_r$S) * zTA$zTA_r)[censor[,4]])/((S_l$S-S_r$S)[censor[,4]])
  
  delta = any(censor[,2])
  
  Q[censor[,2],(p+q+1):(p+q+m)] = (as.numeric(1/h0[censor[,2],]) * psi[censor[,2],]) - (eXBeta[censor[,2]])* (PsiStar_r$Psi_star[censor[,2],])
  Q[censor[,1],(p+q+1):(p+q+m)] = - (eXBeta[censor[,1]])* (PsiStar_r$Psi_star[censor[,1],])
  Q[censor[,3],(p+q+1):(p+q+m)] = (S_r$S[censor[,3]] * eXBeta[censor[,3]] / (1 - S_r$S[censor[,3]])) * PsiStar_r$Psi_star[censor[,3],]
  Q[censor[,4],(p+q+1):(p+q+m)] = (S_r$S[censor[,4]] * eXBeta[censor[,4]] / (S_l$S[censor[,4]] - S_r$S[censor[,4]])) * PsiStar_r$Psi_star[censor[,4],] - 
    (S_l$S[censor[,4]] * eXBeta[censor[,4]] / (S_l$S[censor[,4]] - S_r$S[censor[,4]])) * PsiStar_l$Psi_star[censor[,4],]
  
  n = nrow(dat.baseline)
  Sp = Q-matrix(rep(c(rep(0,p+q),p_t$TwoLRtheta),n),n,byrow=T)/n
  Q = t(Sp)%*%Sp
  
  
  
  M_2 = Hess + 2*lambda*Rstar #m2
  
  Minv_2=corr=matrix(0,p+q+m,p+q+m)
  diag(corr)=rep(1,m+q+p)
  corr[!pos,]=0
  Minv_2[pos,pos]=solve(M_2[pos,pos])
  
  A_eta = corr %*% Minv_2 %*% t(corr)
  
  cov_H = A_eta %*% (Hess) %*% A_eta
  cov_Q = A_eta %*% Q %*% A_eta
  #cov_H=corr%*%(Minv_2%*%Hess%*%Minv_2)%*%t(corr)
  se_H=sqrt(diag(cov_H))
  se_Q=sqrt(diag(cov_Q))
  
  
  
  pAIC = -2*log_lik + 2*p*q*m
  
  out = list(parameters = parameters, theta_info = theta_info, kn = kn, func_est = func_est, iterations = c(iter, it),
             lambda = lambda, conv_record = c(reg.conv, sm.conv), se_H = se_H, 
             covar_H = cov_H, Hess = Hess, grad = grad, covar_Q = cov_Q, se_Q = se_Q,
             pll = log_lik, pAIC = pAIC)
  
  return(out)
  
  
}





