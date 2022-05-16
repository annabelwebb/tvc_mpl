save = matrix(0, nrow = 500, ncol = 26)
save[,1] = c(1:500)

colnames(save) = c("id", "e_p", "r_p", "mpl_b1", "mpl_b2", "mpl_seb1", "mpl_seb2", "mpl_g1", "mpl_seg1",
                   "pl_b1", "pl_b2", "pl_seb1", "pl_seb2", "pl_g1", "pl_seg1", "true_S0_1", 
                   "true_S0_2", "true_S0_3", "mpl_S0_1", "mpl_S0_2", "mpl_S0_3", "pl_S0_1", 
                   "pl_S0_2", "pl_S0_3", "lambda", "mpl_converge")

for(s in 1:500){
  
  #set up data
  sample1 = gen_tvc_right_gomp(1000, c(1,-0.5), c(-1), 0.2, 0.5, 0.6, 1)
  
  
  dat = as.data.frame(sample1$dat)
  dat.baseline = as.data.frame(sample1$dat.baseline)
  dat.left = as.data.frame(sample1$dat.left)
  
  surv.obj = Surv(time = dat$start, time2 = dat$end, type = "counting", event = dat$last_record)
  
  #save censoring proportions
  save[s,2] = sum(dat.baseline$delta == 1)/1000
  save[s,3] = sum(dat.baseline$delta == 0)/1000
  
  ##MPL
  #fit model
  ctrl = tvc_mpl_control(lambda = 0, iter = c(10,500), n_knots = 6, par_initial = c(0,0,1), 
                         range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
  
  mpl.fit = tvc_fit(dat, dat.baseline, dat.left, ctrl)
  
  #save beta info
  save[s,4] = mpl.fit$parameters$beta[1]
  save[s,5] = mpl.fit$parameters$beta[2]
  save[s,6] = mpl.fit$se_H[1]
  save[s,7] = mpl.fit$se_H[2]
  
  #save gamma info
  save[s,8] = mpl.fit$parameters$gamma[1]
  save[s,9] = mpl.fit$se_H[3]
  
  ##PL COX
  #fit model
  pl.fit = coxph(surv.obj ~ dat$x1_long + dat$x2_long + dat$z1, data = dat, id = dat$i_long )
  
  #save beta info
  save[s,10] = pl.fit$coefficients[1]
  save[s,11] = pl.fit$coefficients[2]
  save[s,12] = sqrt(pl.fit$var[1,1])
  save[s,13] = sqrt(pl.fit$var[2,2])
  
  #save gamma info
  save[s,14] = pl.fit$coefficients[3]
  save[s,15] = sqrt(pl.fit$var[3,3])
  
  
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

  alpha = 0.2
  lambda = 0.5
  true_S0 = exp(-((lambda/alpha)*exp(alpha*times.true) -(lambda/alpha)))
  save[s,16:18] = true_S0
  
  #mpl S0
  mpl_S0 = mpl.fit$func_est$S0_est[time.true.ind]
  save[s, 19:21] = mpl_S0
  
  #pl S0
  #pick up times close to quantiles
  pl.surv.fit = survfit(pl.fit)
  times.true.pl = NULL
  time.true.ind.pl = NULL
  for(t in 1:3){
    time.true.ind.pl = c(time.true.ind.pl,max(which(pl.surv.fit$time < time.quant[t])))
    times.true.pl = c(times.true.pl, pl.surv.fit$time[time.true.ind.pl[t]])
  }
  pl_S0 = pl.surv.fit$surv[time.true.ind.pl]
  save[s,22:24] = pl_S0
  
  
  #save extra info on mpl
  save[s,25] = mpl.fit$lambda
  save[s,26] = mpl.fit$conv_record[1]
  
  
  print(c(s,save[s,26], save[s,3], mean(save[1:s, 8]), mean(save[1:s, 9]), sd(save[1:s, 8]), mean(save[1:s, 14])) )
  
  
  
}


write.csv(save,"right_censoring/04_n1000_c075_6k.csv")



