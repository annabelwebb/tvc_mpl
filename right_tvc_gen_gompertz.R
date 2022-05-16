gen_tvc_right_gomp = function(n, b, g, tau1, tau2, c1, c2){
  id = c(1:n)
  i_long = sort(c(id,id))
  
  
  #generate time fixed covariates X
  x1=rbinom(n,1,0.5)
  x2=runif(n)
  X = as.matrix(cbind(x1,x2), ncol = 2)
  
  #generate series of time varying covariates Z
  z1 = rep(c(0,1),n)
  Z = as.matrix(z1, ncol = 1)
  
  #generate a change point time for each i
  change = runif(n, tau1, tau2)
  
  #generate event times
  U_Y = runif(n)
  log.u=log(U_Y)
  beta=matrix(b,ncol = 1)
  gamma=matrix(g, ncol = 1)
  y_i=rep(NA,n)
  alpha = 0.2
  lambda = 0.5
  
  for(i in 1:n){
    H_i = (lambda * exp(X[i,]%*%beta) / alpha) * (exp(alpha*change[i]) - 1)
    
    if( (-log.u[i]) < H_i){
      y_i[i] = (1/alpha) * log(1 + alpha* (-log.u[i]) / (lambda * exp(X[i,]%*%beta))  )
    }else{
      y_i[i] = (1/alpha) * log( (alpha*(-log.u[i])/(lambda * exp(X[i,]%*%beta + gamma)))  - ((exp(alpha*change[i]) - 1 - exp(gamma + alpha*change[i]))/ (exp(gamma)))   )
    }
  }
  
  #generate censoring times
  c_i = runif(n, c1, c2)
  delta = as.numeric(y_i < c_i)
  TL = y_i
  TR = y_i
  TL[delta == 0] = c_i[delta == 0]
  TR[delta == 0] = Inf
  
  #create long format data
  TL_long = TR_long = x1_long = x2_long = delta_long = rep(0, length(i_long))
  
  
  
  for(i in 1:n){
    ind = which(i_long == i)
    x1_long[ind] = x1[i]
    x2_long[ind] = x2[i]
    TL_long[ind] = TL[i]
    TR_long[ind] = TR[i]
    delta_long[ind] = delta[i]
  }
  
  start = rep(0, length(i_long))
  end = rep(0, length(i_long))
  last_record = rep(0, length(i_long))
  
  dat = cbind(i_long, x1_long, x2_long, z1, start, end, TL_long, TR_long, delta_long, last_record)
  
  #create start and end columns
  for(i in 1:n){
    ind = which(i_long == i)
    if(TL[i] < change[i]){ #no change in Z observed
      dat[ind[1],6] = TL[i]
      dat[ind[2],5] = NA
    }else{ #change in Z observed
      dat[ind[1],6] = change[i]
      dat[ind[2],5] = change[i]
      dat[ind[2],6] = TL[i]
    }
  }
  
  if(any(is.na(dat[,5]))){
    dat = dat[-which(is.na(dat[,5])),]
  }
  
  #identify last record
  
  
  for(i in 1:n){
    ind = which(dat[,1] == i)
    max.ind = max(ind)
    dat[max.ind,10] = 1
  }
  
  
  
  rep = (as.data.frame(dat) %>% group_by(i_long) %>% tally())$n
  
  i = id
  
  dat.left = dat
  dat.left[,-1]=0
  
  dat.baseline = cbind(i, x1, x2, TL, TR, delta, rep)
  out = list(dat = dat, dat.baseline = dat.baseline, dat.left=dat.left)
  
  return(out)
  
}







