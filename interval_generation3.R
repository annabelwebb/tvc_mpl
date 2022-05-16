gen_tvc_interval = function(n, b, g, pi_E, t1, t2, a1, a2){
  #n = sample size
  #b = true beta value (list of length 2)
  #g = true gamma value (list of length 1)
  #t1, t2 = time interval within which z will change
  #a1, a2 = interval censoring width adjustment
  
  #create id columns
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
  change = runif(n, t1, t2)
  
  #generate event times from Cox model
  #h0(t) = 3t^2 and H0(t) = t^3
  U_Y = runif(n)
  neg.log.u=-log(U_Y)
  beta=matrix(b,ncol = 1)
  gamma=matrix(g, ncol = 1)
  y_i=rep(NA,n)
  
  for(i in 1:n){
    H_i = exp(X[i,]%*%beta)*change[i]^3
    
    if(neg.log.u[i] <= H_i){
      y_i[i] = (neg.log.u[i]/(exp(X[i,]%*%beta)))^(1/3)
    }else{
      y_i[i] = ((neg.log.u[i] - exp(X[i,]%*%beta)*change[i]^3 + exp(X[i,]%*%beta + gamma)*change[i]^3)/exp(X[i,]%*%beta + gamma))^(1/3)
    }
  }
  
  #generate censoring type, (TL, TR), and start, end columns
  TL=TR=delta=rep(NA,n)
  
  start=end=last_record=left_position = tau = rep(0,n*2)
  #start and end columns give times that z changes
  #last_record is an indicator of value 1 if that row is the last row for a given i, 0 otherwise
  #left_position is an indicator of value 1 if the (start, end) interval for that row contains TL (for interval censored i)
  #tau is an indicator of value 1 if the "start" time for that row is > TL for that row (for interval censored i)
  
  for(i in 1:n){
    ind = which(i_long == i)
    
    U_E = runif(1)
    U_L = runif(1)
    U_R = runif(1)
    
    if(U_E < pi_E){ #event time
      TL[i]=y_i[i]
      TR[i]=y_i[i]
      delta[i] = 1
      
      #create start and end columns
      if(change[i] < y_i[i]){ #two values for z
        end[ind[1]] = change[i]
        start[ind[2]] = change[i]
        end[ind[2]] =  y_i[i]
        last_record[ind[2]] = 1
      }else{ #one value for z
        end[ind[1]] = TL[i]
        end[ind[2]] = NA
        last_record[ind[1]] = 1
      }
      
    }else{
      if(a1*U_L <= y_i[i] & y_i[i] <= a2*U_R){ #interval censoring
        TL[i] = a1*U_L
        TR[i] = a2*U_R
        delta[i] = 3
        
        #create start and end columns
        if(change[i] < TR[i]){ #two values of z
          end[ind[1]] = change[i]
          start[ind[2]] = change[i]
          end[ind[2]] = TR[i]
          last_record[ind[2]] = 1
        }else{
          end[ind[1]] = TR[i]
          end[ind[2]] = NA
          last_record[ind[1]] = 1
        }
        
        
      }else if(a2*U_R < y_i[i]){ #right censoring
        TL[i] = a2*U_R
        TR[i] = Inf
        delta[i] = 0
        
        #create start and end columns
        if(change[i]<TL[i]){ #two values for z
          end[ind[1]] = change[i]
          start[ind[2]] = change[i]
          end[ind[2]] = TL[i]
          last_record[ind[2]] = 1
        }else{ #one values for z
          end[ind[1]] = TL[i]
          end[ind[2]] = NA
          last_record[ind[1]] = 1
        }
        
      }else if(y_i[i] < a1*U_L){ #left censoring
        TL[i] = -Inf
        TR[i] = a1*U_L
        delta[i] = 2
        
        #create start and end columns
        if(change[i]<TR[i]){ #two values for z
          end[ind[1]] = change[i]
          start[ind[2]] = change[i]
          end[ind[2]] = TR[i]
          last_record[ind[2]] = 1
        }else{ #one value for z
          end[ind[1]] = TR[i]
          end[ind[2]] = NA
          last_record[ind[1]] = 1
        }
        
      }
    }
  }
  
  
  #create long format columns
  TL_long = TR_long = x1_long = x2_long = delta_long = rep(0, length(i_long))
  
  for(i in 1:n){
    ind = which(i_long == i)
    x1_long[ind] = x1[i]
    x2_long[ind] = x2[i]
    TL_long[ind] = TL[i]
    TR_long[ind] = TR[i]
    delta_long[ind] = delta[i]
  }
  
  #long format data frame
  dat = cbind(i_long, x1_long, x2_long, z1, start, end, TL_long, TR_long, delta_long, last_record, left_position, tau)
  
  #delete all irrelevant rows from long format data
  #(irrelevant = the second value of z was not actually observed because event time/censoring time was before the change point time)
  if(any(is.na(dat[,6]))){
    dat = dat[-which(is.na(dat[,6])),]
  }
  
  #assign value of 1 to left_position or tau columns where relevant
  for(i in 1:n){
    if(delta[i] == 3){
      ind = which(dat[,1] == i)
      
      if(TL[i] < dat[ind[1],6]){
        dat[ind[1],11] = 1
        dat[ind[1],12] = 1
      }else{
        dat[ind[2],11] = 1
        dat[ind[1],12] = 1
        dat[ind[2],12] = 1
        
      }
    }
  }
  
  #create long format data that ends at TL for interval censored i's
  dat.left = dat
  dat.left[which(dat.left[,12]!=1),-1] = 0
  
  left.i = which(delta==3)
  
  for(li in left.i){
    ind = which(dat.left[,1] == li)
    if(dat.left[ind[1],11] == 1){ #left time in first interval
      dat.left[ind[1],6] = dat.left[ind[1],7]
      dat.left[ind[2],-1] = 0
      
    }else{ #left time in second interval
      dat.left[ind[2],6] = dat.left[ind[2],7]
    }
  }
  
  #how many observations for each i
  rep = (as.data.frame(dat) %>% group_by(i_long) %>% tally())$n
  
  #create baseline dataset
  i = id
  dat.baseline = cbind(i, x1, x2, TL, TR, delta, rep)
  
  out = list(dat = dat, dat.baseline = dat.baseline, dat.left=dat.left)
  
  
  return(out)
  
}

