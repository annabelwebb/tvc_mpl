
gen_int_compare = function(n, b, g, t1, t2){
  
  #create id column
  id_long = sort(rep(c(1:n), 4))
  
  #create covariates
  x1=rbinom(n,1,0.5)
  x2=runif(n)
  
  X = as.matrix(cbind(x1,x2), ncol = 2)
  z1 = rep(c(0,1,0,1),n)
  
  #generate random changing times for each i
  exam_times = rep(0, length(id_long))
  for(i in 1:n){
    ind.long = which(id_long == i)
    change_i = cumsum(runif(3, t1, t2))
    exam_times[ind.long[1:3]] = change_i
    
  }
  
  #generate event times for each i
  U_Y = runif(n)
  neg.log.u=-log(U_Y)
  beta=matrix(b,ncol = 1)
  gamma=matrix(g, ncol = 1)
  y_i=rep(NA,n)
  
  for(i in 1:n){
    ind.long = which(id_long == i)
    change_i = exam_times[ind.long[1:3]]
    
    H_t1 = exp(X[i,]%*%beta) * change_i[1]^3
    H_t2 = exp(X[i,]%*%beta) * change_i[1]^3 + exp(X[i,]%*%beta + gamma)*change_i[2]^3 - exp(X[i,]%*%beta + gamma)*change_i[1]^3
    H_t3 = exp(X[i,]%*%beta)*change_i[1]^3 + exp(X[i,]%*%beta + gamma)*change_i[2]^3 - exp(X[i,]%*%beta + gamma)*change_i[1]^3 +
      exp(X[i,]%*%beta)*change_i[3]^3 - exp(X[i,]%*%beta)*change_i[2]^3
    
    if(neg.log.u[i] < H_t1){
      y_i[i] = (neg.log.u[i]/exp(X[i,]%*%beta))^(1/3)
    }else if(H_t1 <= neg.log.u[i] & neg.log.u[i] < H_t2){
      y_i[i] = ((neg.log.u[i] - exp(X[i,]%*%beta)* (change_i[1]^3 - exp(gamma)*change_i[1]^3))/exp(X[i,]%*%beta + gamma) )^(1/3)
    }else if(H_t2 <= neg.log.u[i] & neg.log.u[i] < H_t3){
      y_i[i] = ((neg.log.u[i] + exp(X[i,]%*%beta)*(1 - exp(gamma)*(change_i[2]^3 - change_i[1]^3)))/exp(X[i,]%*%beta))^(1/3)
    }else if(H_t3 <= neg.log.u[i]){
      y_i[i] = ((neg.log.u[i] - exp(X[i,]%*%beta)*(1 - exp(gamma))*(change_i[1]^3 - change_i[2]^3 + change_i[3]^3) )/ exp(X[i,]%*%beta + gamma) )^(1/3)
      
    }
  }
  
  
  #figure out where each y_i fits into the exam times and create status column
  id_long_final = status = times = z.obs = NULL
  for(i in 1:n){
    ind.long = which(id_long == i)
    
    assessment.n = ceiling(runif(1,1,4))
    assessment.times = cumsum(runif(assessment.n, 0.1,0.3))
    
    if(max(assessment.times)< y_i[i]){ #right censored
      status = c(status, rep(0,length(assessment.times)))
      id_long_final = c(id_long_final, rep(i, length(assessment.times)))
      times = c(times, assessment.times)
      
      change_i = exam_times[ind.long[1:3]]
      z_i = rep(0, assessment.n)
      z_i[which(assessment.times < change_i[1])] = 0
      z_i[which(change_i[1] <= assessment.times & assessment.times < change_i[2])] = 1
      z_i[which(change_i[2] <= assessment.times & assessment.times < change_i[3])] = 0
      z_i[which(assessment.times >= change_i[3])] = 1
      
      z.obs = c(z.obs, z_i)
      
    }else{
      last.assessment = min(assessment.times[which(assessment.times>y_i[i])])
      last.as.ind = which(assessment.times == last.assessment)
      
      status = c(status, rep(0, (last.as.ind -1)),1)
      id_long_final = c(id_long_final, rep(i,last.as.ind))
      times = c(times, assessment.times[1:last.as.ind])
      
      change_i = exam_times[ind.long[1:3]]
      z_i = rep(0, length(assessment.times[1:last.as.ind]))
      z_i[which(assessment.times[1:last.as.ind] < change_i[1])] = 0
      z_i[which(change_i[1] <= assessment.times[1:last.as.ind] & assessment.times[1:last.as.ind] < change_i[2])] = 1
      z_i[which(change_i[2] <= assessment.times[1:last.as.ind] & assessment.times[1:last.as.ind] < change_i[3])] = 0
      z_i[which(assessment.times[1:last.as.ind] >= change_i[3])] = 1
      
      z.obs = c(z.obs, z_i)
      
    }
    
    
  }
  
  x1_long = x2_long = rep(0, length(id_long_final))
  for(i in 1:n){
    ind.long = which(id_long_final == i)
    x1_long[ind.long] = x1[i]
    x2_long[ind.long] = x2[i]
  }
  
  
  
  dat = data.frame(cbind(id_long_final, times, z.obs, status, x1_long, x2_long))
  colnames(dat) = c("id", "times", "z", "status", "x1", "x2")
  
  
  
  #adjust above data frame to make it suitable for my model
  
  #create adjusted 'dat' frame
  
  dat.adj = dat
  dat.adj$start = dat.adj$delta_long = dat.adj$left_position = dat.adj$last_record = dat.adj$tau = 
    dat.adj$TL_long = dat.adj$TR_long = rep(0, nrow(dat))
  n = max(dat.adj$id)
  for(i in 1:n){
    ind.long = which(dat.adj$id == i)
    dat.adj$start[ind.long[-1]] = dat.adj$times[ind.long[-length(ind.long)]]
    
    dat.adj$last_record[max(ind.long)] = 1
    
    if(any(dat.adj$status[ind.long] == 1)){
      if(length(ind.long)> 1){
        dat.adj$TL_long[ind.long] = dat.adj$times[ind.long[(length(ind.long)-1)]]
        dat.adj$TR_long[ind.long] = dat.adj$times[max(ind.long)]
        
        dat.adj$delta_long[ind.long] = 3
        dat.adj$left_position[ind.long[(length(ind.long)-1)]] = 1
        dat.adj$tau[ind.long[1:(length(ind.long)-1)]] = 1
      }else{
        dat.adj$TL_long[ind.long] = -Inf
        dat.adj$TR_long[ind.long] = dat.adj$times[max(ind.long)]
        
        dat.adj$delta_long[ind.long] = 2
      }
      
    }else{
      dat.adj$TL_long[ind.long] = dat.adj$times[max(ind.long)]
      dat.adj$TR_long[ind.long] = Inf
    }
    
    
    
  }
  
  
  id = c(1:max(dat.adj$id))
  x1 = x2 = TL = TR = delta = rep(0, n)
  for(i in 1:n){
    ind.long = which(dat.adj$id == i)
    
    x1[i] = dat.adj$x1[ind.long[1]]
    x2[i] = dat.adj$x2[ind.long[1]]
    TL[i] = dat.adj$TL_long[ind.long[1]]
    TR[i] = dat.adj$TR_long[ind.long[1]]
    delta[i] = dat.adj$delta_long[ind.long[1]]
    
  }
  
  dat.baseline.adj = data.frame(cbind(id, x1, x2, TL, TR, delta))
  dat.baseline.adj$rep = (as.data.frame(dat.adj) %>% group_by(id) %>% tally())$n
  
  
  #create long format data that ends at TL for interval censored i's
  dat.left.adj = dat.adj
  dat.left.adj[which(dat.left.adj$tau !=1),-1] = 0
  
  left.i = which(dat.baseline.adj$delta==3)
  
  for(li in left.i){
    #replace "end" time with left interval time
    #make everything after left interval time 0
    ind = which(dat.left.adj[,1] == li)
    left.pos.ind = ind[which(dat.left.adj$left_position[ind] == 1)]
    dat.left.adj$times[left.pos.ind] = dat.left.adj$TL_long[left.pos.ind] 
    #dat.left[left.pos.ind,4] = log(dat.left[left.pos.ind,6])
    
  }
  
  dat.adj2 = data.frame(dat.adj$id, dat.adj$x1, dat.adj$x2, dat.adj$z, dat.adj$start, dat.adj$times, 
                        dat.adj$TL_long, dat.adj$TR_long, dat.adj$delta_long, dat.adj$last_record, 
                        dat.adj$left_position, dat.adj$tau)
  dat.left.adj2 = data.frame(dat.left.adj$id, dat.left.adj$x1, dat.left.adj$x2, dat.left.adj$z, dat.left.adj$start, dat.left.adj$times, 
                             dat.left.adj$TL_long, dat.left.adj$TR_long, dat.left.adj$delta_long, dat.left.adj$last_record, 
                             dat.left.adj$left_position, dat.left.adj$tau)
  
  colnames(dat.adj2) = c("i_long", "x1_long", "x2_long", "z1", "start", "end", "TL_long", "TR_long",
                         "delta_long", "last_record", "left_position", "tau")
  
  
  
  colnames(dat.left.adj2) =  c("i_long", "x1_long", "x2_long", "z1", "start", "end", "TL_long", "TR_long",
                               "delta_long", "last_record", "left_position", "tau")
  
  
  out = list(dat = dat, dat.adj2 = dat.adj2, dat.baseline.adj = dat.baseline.adj,
             dat.left.adj2 = dat.left.adj2)
  
}



