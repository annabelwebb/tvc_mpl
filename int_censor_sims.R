save = matrix(0, nrow = 500, ncol = 22)
save[,1] = c(1:500)

colnames(save) = c("id", "e_p", "r_p", "l_p", "i_p",
                   "mpl_b1", "mpl_b2", "mpl_seb1", "mpl_seb2", "mpl_g1", "mpl_seg1", 
                   "true_S0_1", "true_S0_2", "true_S0_3", "mpl_S0_1", "mpl_S0_2", "mpl_S0_3", 
                    "mpl_S0_1_se" , "mpl_S0_2_se", "mpl_S0_3_se" ,"lambda", "mpl_converge")

for(s in 1:500){
  
  #set up data
  sample1 = gen_tvc_interval(200, c(-1, 0.5), c(0.5), 0.3, 0.5, 1, 0.8, 1.2)
  
  dat = as.data.frame(sample1$dat)
  dat.baseline = as.data.frame(sample1$dat.baseline)
  dat.left = as.data.frame(sample1$dat.left)
  
  #save censoring proportions
  save[s,2] = sum(dat.baseline$delta == 1)/200
  save[s,3] = sum(dat.baseline$delta == 0)/200
  save[s,4] = sum(dat.baseline$delta == 2)/200
  save[s,5] = sum(dat.baseline$delta == 3)/200
  
  ##MPL
  #fit model
  ctrl = tvc_mpl_control(lambda = 0, iter = c(10,1000), n_knots = 2, par_initial = c(0,0,1), 
                         range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
  
  
  
  try=try(tvc_fit(dat, dat.baseline, dat.left, ctrl))
  
  if(class(try)!="try-error"){
    mpl.fit=try
  }else if(class(try)=="try-error"){
    while(class(try)=="try-error"){
      sample1 = gen_tvc_interval(200, c(-1, 0.5), c(0.5), 0.3, 0.5, 1, 0.8, 1.2)
      
      dat = as.data.frame(sample1$dat)
      dat.baseline = as.data.frame(sample1$dat.baseline)
      dat.left = as.data.frame(sample1$dat.left)
      
      #save censoring proportions
      save[s,2] = sum(dat.baseline$delta == 1)/200
      save[s,3] = sum(dat.baseline$delta == 0)/200
      save[s,4] = sum(dat.baseline$delta == 2)/200
      save[s,5] = sum(dat.baseline$delta == 3)/200
      
      ##MPL
      #fit model
      ctrl = tvc_mpl_control(lambda = 0, iter = c(10,1000), n_knots = 2, par_initial = c(0,0,1), 
                             range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
      
      try=try(tvc_fit(dat, dat.baseline, dat.left, ctrl))
    }
    mpl.fit=try
  }
  
  
  #save beta info
  save[s,6] = mpl.fit$parameters$beta[1]
  save[s,7] = mpl.fit$parameters$beta[2]
  save[s,8] = mpl.fit$se_H[1]
  save[s,9] = mpl.fit$se_H[2]
  
  #save gamma info
  save[s,10] = mpl.fit$parameters$gamma[1]
  save[s,11] = mpl.fit$se_H[3]
  
  
  ##BASELINE SURVIVAL ESTIMATION
  #time and true function
  time.quant = quantile(mpl.fit$func_est$v)[2:4]
  #pick up real times closest to these quantiles
  times.true = NULL
  time.true.ind = NULL
  for(t in 1:3){
    time.true.ind = c(time.true.ind,max(which(mpl.fit$func_est$v < time.quant[t])))
    times.true = c(times.true, mpl.fit$func_est$v[time.true.ind[t]])
  }
  
  true_S0 = exp(-times.true^3)
  save[s,12:14] = true_S0
  
  #mpl S0
  mpl_S0 = mpl.fit$func_est$S0_est[time.true.ind]
  save[s, 15:17] = mpl_S0
  
  #save asymptotic SE
  cov.theta = mpl.fit$covar_H[-c(1:3), -c(1:3)]
  
  Psi.t = Psi_f(events = c(times.true[1]), knts = mpl.fit$kn$int_knots, Boundary.knots = mpl.fit$kn$bound_knots)
  S.t = (mpl.fit$func_est$S0_est[time.true.ind])[1]
  save[s, 18] = sqrt(t(t(Psi.t) %*% S.t) %*% cov.theta %*% (t(Psi.t) %*% S.t))
  
  Psi.t = Psi_f(events = c(times.true[2]), knts = mpl.fit$kn$int_knots, Boundary.knots = mpl.fit$kn$bound_knots)
  S.t = (mpl.fit$func_est$S0_est[time.true.ind])[2]
  save[s, 19] = sqrt(t(t(Psi.t) %*% S.t) %*% cov.theta %*% (t(Psi.t) %*% S.t))
  
  Psi.t = Psi_f(events = c(times.true[3]), knts = mpl.fit$kn$int_knots, Boundary.knots = mpl.fit$kn$bound_knots)
  S.t = (mpl.fit$func_est$S0_est[time.true.ind])[3]
  save[s, 20] = sqrt(t(t(Psi.t) %*% S.t) %*% cov.theta %*% (t(Psi.t) %*% S.t))
  
  
  #save extra info on mpl
  save[s,21] = mpl.fit$lambda
  save[s,22] = mpl.fit$conv_record[1]
  
  
  print(c(s,save[s,22], save[s,2], mean(save[1:s, 10]), mean(save[1:s, 11]), sd(save[1:s, 10]), 
          mean(save[1:s, 12] - save[1:s, 15]), mean(save[1:s, 18]), sd(save[1:s, 15]) ))
  
  
  
}



write.csv(save,"int_censoring/02_n200_c075_2k_SE_adjusted.csv")
   



