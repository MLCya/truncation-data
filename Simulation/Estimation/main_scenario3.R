
#####    main scenrio 3:  F&G exponential distribution
###    case1： z ~ U[0,5]:  beta <- c(0.3,1) ;                       rate.G <- 10;  OP_1 =  0.7914179
###    case2： z ~ U[0,5]:  beta <- c(0.1, 0.9, 0.07);               rate.G <- 7;   OP_2 =  0.7286483
###    case3： z1 ~ U[0,5], z2 ~ B(1,0.5):  beta <- c(0.1,1.2, 0.1); rate.G <- 10;  OP_3 =  0.7740895
##### modelling:  theta  = exp(beta %*% qz)


library(rootSolve)
library(cubature)
eeps <-  1.490116e-08  



digam <- function(n_big, n_small){ 
  out=0
  if(n_small>0 &  n_big>=n_small)  out=sum(log(n_big-n_small+c(1:n_small)))    
  out
}


#####  define a new log function  log_n    ###
log_n <- function(z,n){
  an=n^2
  ind=(z<1/an)
  y1 =log(1/an)-1.5+2*an*z-(an*z)^2/2
  y1[!ind]=0
  z[ind] = 1
  y2=log(z)
  y=y1+y2 
  return(y)
}

###### given alpha to estimate N
MaxN.alpha <- function(n_small, alpha){  
  myfun <- function(N){  
    tmp = N +1 - (1:n_small)
    sum(1/tmp)+log(alpha)  
  }
  multiroot(f = myfun,start = c(n_small),maxiter = 5000)$root
  #lower =  n_small
  #upper =  -n_small/log(alpha)  + n_small - 1
  #uniroot(myfun, c(lower+2, upper),extendInt = "yes",maxiter = 50000)$root ## ,extendInt = "yes"
}



f_b <- function(x,theta){
  1 - exp(-theta*x)
}

f_s <- function(x,theta){
  theta*exp(-theta*x) 
}

log_fs <- function(x,theta){
  log(theta) - theta*x
}

log_fb <- function(x,theta){## log(1-F(y,theta))
  -theta*x
}

f1_s <- function(x,theta){
  (1-theta*x)*exp(-theta*x)
}

f2_s <- function(x,theta){
  (theta*x - 2)*x*exp(-theta*x)
}


f1_b <- function(x,theta){
  x*exp(-theta*x)
}

f2_b <- function(x,theta){
  -(x^2)*exp(-theta*x)
}






##---------  Calculate the sampling probability  ------------###
###   case1 :  z ~ U[0,10] ,true theta : h(z,beta) = beta0 +  beta1*z
O_P1 <- function(beta,rate.G){  
  comf <- function(yz){#yz  = (y,z) 
    theta <- beta[1] + beta[2]*yz[2]
    f_b(yz[1],theta = theta)*f_s(yz[1],rate.G)*dunif(yz[2],min = 0, max = 5)
  }
  1 - adaptIntegrate(comf,lowerLimit = c(0,0),upperLimit = c(Inf,5),tol = 1e-05)$integral
}

###  case 2 :  z ~ U[0,5], true theta : h(z,beta) = beta0 +  beta1*z + beta2*z^2
O_P2  <- function(beta,rate.G){  
  comf <- function(yz){#yz  = (y,z) 
    theta <- beta[1] + beta[2]*yz[2] + beta[3]*yz[2]^2
    f_b(yz[1],theta = theta)*f_s(yz[1],rate.G)*dunif(yz[2],min = 0,max = 5)
  }
  1 - adaptIntegrate(comf,lowerLimit = c(0,0),upperLimit = c(Inf,5),tol = 1e-05)$integral
}



### case3 : z1 ~ U[0,10], z2 ~ B(1,0.5), true theta : h(z,beta) = beta0 +  beta1*z + beta2*z2
O_P3 <- function(beta,rate.G){
  comf <- function(yz){
    y <- yz[1]
    z1 <- yz[2]
    z2 <- yz[3]
    theta0 <- beta[1]+beta[2]*z1        #z2=0
    theta1 <- beta[1]+beta[2]*z1+beta[3]#z2=1
    
    (f_b(y, theta0) + f_b(y,theta1) )*f_s(y,rate.G)/5/2
  }
  1 - adaptIntegrate(comf,lowerLimit = c(0,0),upperLimit = c(Inf,5))$integral
}






##-------------          generating  data: -----------#
##  case1:
gen1_obs <- function(beta,rate.G,n.big){#beta is a vector of 2 dimension,qz=(1,z)
  y = rexp(n.big,rate = rate.G)
  z0 = rep(1,n.big)
  z = runif(n.big,min = 0,max = 5)
  theta_t = beta[1]*z0+beta[2]*z
  x = rexp(n.big,rate = theta_t)
  
  ind = (x > y)
  obs <- data.frame(obs.x= x[ind],obs.y=y[ind],obs.z0=z0[ind],obs.z=z[ind],theta.obs = theta_t[ind])
  return(obs)
  
}



##  case2 :
gen2_obs <- function(beta,rate.G,n.big){#beta is a vector of 3 dimension,qz=(1,z,z^2)
  y = rexp(n.big,rate = rate.G)
  z0 = rep(1,n.big)
  z  = runif(n.big,min = 0,max = 5)
  theta_t = beta[1]*z0+beta[2]*z + beta[3]*z^2
  x = rexp(n.big,rate = theta_t)
  
  ind = (x > y)
  obs <- data.frame(obs.x= x[ind],obs.y=y[ind],obs.z0=z0[ind],obs.z=z[ind],theta.obs = theta_t[ind])
  return(obs)
}

##  caes 3 :
gen3_obs <- function(beta,rate.G, n.big){
  y <- rexp(n.big,rate.G)
  z0 <- rep(1,n.big)
  z1 <- runif(n.big,0,5)
  z2  = NULL
  for (i in  1:n.big) {
    z2[i] = rbinom(1,1,0.5)
  } 
  
  theta_t <- beta[1] + beta[2]*z1 + beta[3]*z2
  x <- NULL
  for (i in 1:n.big) {
    x[i] <- rexp(1,rate = theta_t[i])
  }
  ind <- ( x > y)
  obs <- data.frame(x.obs = x[ind],y.obs = y[ind],z0.obs = z0[ind],z1.obs = z1[ind],z2.obs = z2[ind],theta.obs = theta_t[ind])
  return(obs)
}





##--------------------    full likelihood procedure     -----------------------#
##    case 1  
maxlike1.par <- function(obs){
  n <- nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  
  
  prof.lik <- function(beta){  ###profile loglikelihood of beta
    theta   <- exp(beta[1]*z0 + beta[2]*z) 
    fbig    <- f_b(y,theta)
    fsmall  <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta) ## log(f(x,theta))
    
    
    fn.iner <- function(gam){ ## gam is a reparameter
      ## reparametrics
      alpha <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      N  <- MaxN.alpha(n_small = n,alpha = alpha )
      
      d <- fbig - alpha
      
      lab = 0  
      if(max(d) * min(d) < 0){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      
      fn <- digam(N,n) + (N-n)*log(alpha) - sum(log(1+lab*d))
      return(-fn )
    }
    optimize(fn.iner,c(0,1))$objective - sum(logfs)
  }
  
  ## the eatimator of beta 
  out <- optim(c(1,1), fn=prof.lik)
  beta_Phat <-  out$par
  
  
  
  theta_Phat  <- exp(beta_Phat[1] + beta_Phat[2]*z )  ##the eatimator of theta = h(z,beta)
  fbig_Phat   <- f_b(y,theta_Phat)
  fsmall_Phat <- f_s(x,theta_Phat)
  
  fbig.1_Phat <- min(fbig_Phat)
  fbig.n_Phat <- max(fbig_Phat)
  logfs_Phat <- log_fs(x,theta_Phat)
  
  fn.outer <- function(gam){ ## ng = c(n.big,gama) is parametric vector
    alpha_ter <- 1/(1+gam^2) * fbig.1_Phat + gam^2/(1+gam^2)*fbig.n_Phat
    N <- MaxN.alpha(n_small = n, alpha = alpha_ter)
    d <- fbig_Phat  - alpha_ter
    
    lab = 0
    if(max(d)*min(d)<0){
      fl  <- function(lab){ sum(d/(1+lab*d))   }  
      low= -1/max(d) + eeps
      up = - 1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- digam(N,n) + (N-n)*log(alpha_ter) - sum(log(1+lab*d))
    return(-f)
  }
  
  result <- optimize(fn.outer,interval = c(-1,1))
  gama_Phat   <- result$minimum
  alpha_Phat <- 1/(1+gama_Phat^2)*fbig.1_Phat + gama_Phat^2/(1+gama_Phat^2)*fbig.n_Phat  ## the eatimator of alpha
  
  N_Phat  <- MaxN.alpha(n,alpha_Phat)       ## the eatimator of N
  #theta_Phat.m <- mean(theta_Phat)
  par_Phat  <- c(N_Phat, alpha_Phat,beta_Phat)
  maxlike   =  -out$value 
  list(par_Phat,maxlike)
}



##-------------    calculating likelihood ratio for parameters    ------------------------------##
### case1 :N=N0 is known 

maxlike1.N <- function(obs,n.big){
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  n <- nrow(obs)
  N <- n.big
  
  prof.lik <- function(beta){#  profile loglikelihood of beta
    theta <- exp(beta[1]*z0 + beta[2]*z)
    fbig <- f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta)
    
    fn.inner <- function(gam){## gam reparametric parameter,gam is scalar
      alpha <- 1/(1+gam^2)*fbig.1 + gam^2/(1+gam^2)*fbig.n
      d <- fbig - alpha
      
      lab = 0
      if(min(d)*max(d) < 0 ){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      fn <- (N - n)*log(alpha) - sum(log(1+lab*d))
      -fn
    }
    optimize(fn.inner,interval = c(-1,1))$objective - sum(logfs )
  }
  maxlike <- digam(N,n) - optim(par = c(1,1),fn=prof.lik)$value  
  maxlike
}



### case1 : alpha = alpha0  is known

maxlike1.A  <- function(obs,alpha){
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z  <- obs[,4]
  n <- nrow(obs)
  n.big = MaxN.alpha(n, alpha) ## the eatimator of N
  
  prof.lik  <- function(beta){  ## profile loglikelihood of beta
    theta <- exp(beta[1]*z0 + beta[2]*z)
    fbig <- f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    logfs <- log_fs(x,theta)
    d = fbig - alpha 
    
    lab = 0
    if(min(d)*max(d) <0 ){
      fl  <- function(lab){ sum(d/(1+lab*d))}  
      low= -1/max(d) + eeps
      up = -1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- sum(logfs) - sum(log( 1 + lab*d))
    -f
  }
  maxlike <- digam(n_big = n.big,n_small = n) + (n.big-n)*log(alpha)- optim(par = c(1,1),fn=prof.lik)$value 
  maxlike
}




##--------------------    conditional likelihood procedure     -----------------------#
##    case 1 
max1.clike <-function(obs){ 
  n = nrow(obs)
  x = obs[, 1]
  y = obs[, 2]
  z0 = obs[,3]
  z = obs[,4]
  
  llc <- function(beta){
    theta <- exp(beta[1]*z0 + beta[2]*z ) 
    logfs <- log_fs(x,theta) 
    logfb <- log_fb(y,theta)
    -sum(logfs - logfb)
  }
  
  beta_tilde = optim(par = c(1,1),fn = llc)$par
  theta_tilde <- exp(beta_tilde[1]*z0+beta_tilde[2]*z)
  
  N_tilde     = sum(1/(1- f_b(y, theta_tilde)))
  alpha_tilde = 1- n/N_tilde
  par_tilde   <- c(N_tilde, alpha_tilde, beta_tilde)
  return(par_tilde)
}




#####################################################################################################
### case 1:  estimate the asymptotical variances of tilde(N), tilde(alpha) and tilde(beta0),tilde(beta1).
Variance1_est <- function(obs){
  x = obs[,1]
  y = obs[,2]
  z0 = obs[,3]
  z = obs[,4]
  qz = cbind(z0,z)
  n = nrow(obs)
  
  til_par = max1.clike(obs)
  til_N = til_par[1]
  til_a = til_par[2]
  til_b0 = til_par[3]
  til_b1 = til_par[4]
  til_theta = exp(til_b0*z0 + til_b1*z)
  
  f1.big <- f1_b(y,theta = til_theta)
  til_phi <- sum(exp( 2*til_theta*y))/til_N
  
  coef_V23 <- NULL 
  coef_V33 <- NULL ## calculate the coefficients of components of V_{33}^{-1}
  for (i in 1:n) {
    fn <- function(x){
      (1/til_theta[i] - 2*x +til_theta[i]*x^2)/exp(til_theta[i]*x)
    }
    int_val <- integrate(fn,lower = y[i],upper = Inf)[[1]]
    sec_val <- f1_b(y[i],til_theta[i])^2*exp(til_theta[i]*y[i])
    coef_V33[i] <- (int_val - sec_val )*til_theta[i]^2*exp(til_theta[i]*y[i])
    
    coef_V23 <- rbind(coef_V23,f1_b(y[i],til_theta[i])*til_theta[i]*exp(2*til_theta[i]*y[i])*qz[i,])
  }  
  til_V23 <-  -apply(coef_V23, 2, sum)/til_N
  V1_est  <- -1/til_N*sum(coef_V33)
  Vz_est  <-  -1/til_N*sum(coef_V33*z)
  Vzs_est <- -1/til_N*sum(coef_V33*z^2)## V_{z^2}
  V_d <- V1_est*Vzs_est - Vz_est^2
  
  ############ the  components of V_{33}^{-1}  ########################
  inV_11 <- Vzs_est/V_d
  inV_12 <- - Vz_est/V_d
  inV_22 <- V1_est/V_d
  com <- c(inV_11,inV_12,inV_12,inV_22)
  
  til_inV33 <- matrix(com,ncol = 2,byrow = T,dimnames = NULL)
  
  Var_N <- til_phi - 1 - sum((til_V23 %*% til_inV33) *  til_V23)
  Var_A <- til_a*(1-til_a)
  Var_b0 <- -inV_11
  Var_b1 <- -inV_22
  
  cbind(Var_N,Var_A,Var_b0,Var_b1)
}






##--------------------    full likelihood procedure     -----------------------#
##    case 2 
maxlike2.par <- function(obs){
  n <- nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  
  
  
  prof.lik  <- function(beta){  ###profile loglikelihood of beta
    theta   <- exp(beta[1]*z0 + beta[2]*z) 
    fbig    <- f_b(y,theta)
    fsmall  <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta) ## log(f(x,theta))
    
    
    fn.iner <- function(gam){ ## gam is a reparameter
      ## reparametrics
      alpha <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      N  <- MaxN.alpha(n_small = n,alpha = alpha )
      
      d <- fbig - alpha
      
      lab = 0  
      if(max(d) * min(d) < 0){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      
      fn <- digam(N,n) + (N-n)*log(alpha) - sum(log(1+lab*d))
      return(-fn )
    }
    optimize(fn.iner,c(0,1))$objective - sum(logfs)
  }
  
  ## the eatimator of beta 
  out <- optim(c(1,1), fn=prof.lik)
  beta_Phat <-  out$par
  
  
  
  theta_Phat  <- exp(beta_Phat[1]*z0 + beta_Phat[2]*z )  ##the eatimator of theta = h(z,beta)
  fbig_Phat   <- f_b(y,theta_Phat)
  fsmall_Phat <- f_s(x,theta_Phat)
  
  fbig.1_Phat <- min(fbig_Phat)
  fbig.n_Phat <- max(fbig_Phat)
  logfs_Phat <- log_fs(x,theta_Phat)
  
  fn.outer <- function(gam){ ## ng = c(n.big,gama) is parametric vector
    alpha_ter <- 1/(1+gam^2) * fbig.1_Phat + gam^2/(1+gam^2)*fbig.n_Phat
    N <- MaxN.alpha(n_small = n, alpha = alpha_ter)
    d <- fbig_Phat  - alpha_ter
    
    lab = 0
    if(max(d)*min(d)<0){
      fl  <- function(lab){ sum(d/(1+lab*d))   }  
      low= -1/max(d) + eeps
      up = - 1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- digam(N,n) + (N-n)*log(alpha_ter) - sum(log(1+lab*d))
    return(-f)
  }
  
  result <- optimize(fn.outer,interval = c(-1,1))
  gama_Phat   <- result$minimum
  alpha_Phat <- 1/(1+gama_Phat^2)*fbig.1_Phat + gama_Phat^2/(1+gama_Phat^2)*fbig.n_Phat  ## the eatimator of alpha
  
  N_Phat  <- MaxN.alpha(n,alpha_Phat)       ## the eatimator of N
  #theta_Phat.m <- mean(theta_Phat)
  par_Phat  <- c(N_Phat, alpha_Phat,beta_Phat)
  maxlike   =  -out$value 
  list(par_Phat,maxlike)
}





##-------------    calculating likelihood ratio for parameters    ------------------------------##
##  case2 :  N=N0 is known 

maxlike2.N <- function(obs,n.big){
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  n <- nrow(obs)
  N <- n.big
  
  prof.lik <- function(beta){#  profile loglikelihood of beta
    theta <- exp(beta[1]*z0 + beta[2]*z)
    fbig <- f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta)
    
    fn.inner <- function(gam){## gam reparametric parameter,gam is scalar
      alpha <- 1/(1+gam^2)*fbig.1 + gam^2/(1+gam^2)*fbig.n
      d <- fbig - alpha
      
      lab = 0
      if(min(d)*max(d) < 0 ){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      fn <- (N - n)*log(alpha) - sum(log(1+lab*d))
      -fn
    }
    optimize(fn.inner,interval = c(-1,1))$objective - sum(logfs )
  }
  maxlike <- digam(N,n) - optim(par = c(1,1),fn=prof.lik)$value  
  maxlike
}



###   case 2: alpha = alpha0  is known

maxlike2.A  <- function(obs,alpha){
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z  <- obs[,4]
  n <- nrow(obs)
  n.big = MaxN.alpha(n, alpha) ## the eatimator of N
  
  prof.lik  <- function(beta){  ## profile loglikelihood of beta
    theta <- exp(beta[1]*z0 + beta[2]*z)
    fbig <- f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    logfs <- log_fs(x,theta)
    d = fbig - alpha 
    
    lab = 0
    if(min(d)*max(d) <0 ){
      fl  <- function(lab){ sum(d/(1+lab*d))}  
      low= -1/max(d) + eeps
      up = -1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- sum(logfs) - sum(log( 1 + lab*d))
    -f
  }
  maxlike <- digam(n_big = n.big,n_small = n) + (n.big-n)*log(alpha)- optim(par = c(1,1),fn=prof.lik)$value 
  maxlike
}






##--------------------    conditional likelihood procedure     -----------------------#
##    case 2

max2.clike <-function(obs){ 
  n = nrow(obs)
  x = obs[, 1]
  y = obs[, 2]
  z0 = obs[,3]
  z = obs[,4]
  
  llc <- function(beta){
    theta <- exp(beta[1]*z0 + beta[2]*z) 
    logfs <- log_fs(x,theta) 
    logfb <- log_fb(y,theta)
    -sum(logfs - logfb)
  }
  
  beta_tilde = optim(par = c(1,1), fn = llc)$par
  theta_tilde <- exp(beta_tilde[1]*z0+beta_tilde[2]*z)
  
  N_tilde     = sum(1/(1- f_b(y, theta_tilde)))
  alpha_tilde = 1- n/N_tilde
  par_tilde   <- c(N_tilde, alpha_tilde, beta_tilde)
  return(par_tilde)
}


#####################################################################################################
### case 2:  estimate the asymptotical variances of tilde(N), tilde(alpha) and tilde(beta0),tilde(beta1).

Variance2_est <- function(obs){
  x = obs[,1]
  y = obs[,2]
  z0 = obs[,3]
  z  = obs[,4]
  qz = cbind(z0,z)
  n = nrow(obs)
  
  til_par = max2.clike(obs)
  til_N = til_par[1]
  til_a = til_par[2]
  til_b0 = til_par[3]
  til_b1 = til_par[4]
  til_theta = exp(til_b0*z0 + til_b1*z)
  
  f1.big <- f1_b(y,theta = til_theta)
  til_phi <- sum(exp( 2*til_theta*y))/til_N
  
  coef_V23 <- NULL 
  coef_V33 <- NULL ## ## calculate the coefficients of components of V_{33}^{-1}
  for (i in 1:n) {
    fn <- function(x){
      (1/til_theta[i] - 2*x +til_theta[i]*x^2)/exp(til_theta[i]*x)
      #f1_s(x, til_theta[i])^2/f_s(x,til_theta[i])
    }
    int_val <- integrate(fn,lower = y[i],upper = Inf)[[1]]
    sec_val <- f1_b(y[i],til_theta[i])^2*exp(til_theta[i]*y[i])
    coef_V33[i] <- (int_val - sec_val)*til_theta[i]^2*exp(til_theta[i]*y[i])
    
    coef_V23 <- rbind(coef_V23,f1_b(y[i],til_theta[i])*til_theta[i]*exp(2*til_theta[i]*y[i])*qz[i,])
  }  
  til_V23 <-  -apply(coef_V23, 2, sum)/til_N
  V1_est  <- -1/til_N*sum(coef_V33)
  Vz_est  <-  -1/til_N*sum(coef_V33*z)
  Vzs_est <- -1/til_N*sum(coef_V33*z^2)## V_{z^2}
  V_d <- V1_est*Vzs_est - Vz_est^2
  
  ############ the  components of V_{33}^{-1}  ########################
  inV_11 <- Vzs_est/V_d
  inV_12 <- - Vz_est/V_d
  inV_22 <- V1_est/V_d
  com <- c(inV_11,inV_12,inV_12,inV_22)
  
  til_inV33 <- matrix(com,ncol = 2,byrow = T,dimnames = NULL)
  
  Var_N <- til_phi - 1 - sum((til_V23 %*% til_inV33) *  til_V23)
  Var_A <- til_a*(1-til_a)
  Var_b0 <- -inV_11
  Var_b1 <- -inV_22
  
  cbind(Var_N,Var_A,Var_b0,Var_b1)
}



##--------------------    full likelihood procedure     -----------------------#
##    case 3 
maxlike3.par <- function(obs){
  n <- nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z1 <- obs[,4]
  z2 <- obs[,5]
  
  
  prof.lik <- function(beta){  ###profile loglikelihood of beta
    theta   <- exp(beta[1]*z0 + beta[2]*z1 + beta[3]*z2) 
    fbig    <- f_b(y,theta)
    fsmall  <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta) ## log(f(x,theta))
    
    
    fn.iner <- function(gam){ ## gam is a reparameter
      ## reparametrics
      alpha <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      N  <- MaxN.alpha(n_small = n,alpha = alpha )
      
      d <- fbig - alpha
      
      lab = 0  
      if(max(d) * min(d) < 0){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      
      fn <- digam(N,n) + (N-n)*log(alpha) - sum(log(1+lab*d))
      return(-fn )
    }
    optimize(fn.iner,c(0,1))$objective - sum(logfs)
  }
  
  ## the eatimator of beta 
  out <- optim(c(1,1,1), fn=prof.lik)
  beta_Phat <-  out$par
  
  
  theta_Phat  <- exp(beta_Phat[1]*z0 + beta_Phat[2]*z1 + beta_Phat[3]*z2)  ##the eatimator of theta = h(z,beta)
  fbig_Phat   <- f_b(y,theta_Phat)
  fsmall_Phat <- f_s(x,theta_Phat)
  
  fbig.1_Phat <- min(fbig_Phat)
  fbig.n_Phat <- max(fbig_Phat)
  logfs_Phat <- log_fs(x,theta_Phat)
  
  fn.outer <- function(gam){ ## ng = c(n.big,gama) is parametric vector
    alpha_ter <- 1/(1+gam^2) * fbig.1_Phat + gam^2/(1+gam^2)*fbig.n_Phat
    N <- MaxN.alpha(n_small = n, alpha = alpha_ter)
    d <- fbig_Phat  - alpha_ter
    
    lab = 0
    if(max(d)*min(d)<0){
      fl  <- function(lab){ sum(d/(1+lab*d))   }  
      low= -1/max(d) + eeps
      up = - 1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- digam(N,n) + (N-n)*log(alpha_ter) - sum(log(1+lab*d))
    return(-f)
  }
  
  result <- optimize(fn.outer, interval = c(-1,1))
  gama_Phat   <- result$minimum
  alpha_Phat <- 1/(1+gama_Phat^2)*fbig.1_Phat + gama_Phat^2/(1+gama_Phat^2)*fbig.n_Phat  ## the eatimator of alpha
  
  N_Phat  <- MaxN.alpha(n,alpha_Phat)       ## the eatimator of N
  #theta_Phat.m <- mean(theta_Phat)
  par_Phat  <- c(N_Phat, alpha_Phat,beta_Phat)
  maxlike   =  -out$value 
  list(par_Phat,maxlike)
}




##-------------    calculating likelihood ratio for parameters    ------------------------------##
###  case3 :  N=N0 is known 

maxlike3.N <- function(obs,n.big){
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z1 <- obs[,4]
  z2 <- obs[,5]
  n <- nrow(obs)
  N <- n.big
  
  prof.lik <- function(beta){#  profile loglikelihood of beta
    theta <- exp(beta[1]*z0 + beta[2]*z1 + beta[3]*z2)
    fbig <- f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta)
    
    fn.inner <- function(gam){## gam reparametric parameter,gam is scalar
      alpha <- 1/(1+gam^2)*fbig.1 + gam^2/(1+gam^2)*fbig.n
      d <- fbig - alpha
      
      lab = 0
      if(min(d)*max(d) < 0 ){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      fn <- (N - n)*log(alpha) - sum(log(1+lab*d))
      -fn
    }
    optimize(fn.inner,interval = c(-1,1))$objective - sum(logfs )
  }
  maxlike <- digam(N,n) - optim(par = c(1,1,1),fn=prof.lik)$value  
  maxlike
}



###   case3 : alpha = alpha0  is known

maxlike3.A  <- function(obs,alpha){
  x <- obs[,1]
  y <- obs[,2]
  z0  <- obs[,3]
  z1  <- obs[,4]
  z2  <- obs[,5]
  n <- nrow(obs)
  n.big = MaxN.alpha(n, alpha) ## the eatimator of N
  
  prof.lik  <- function(beta){  ## profile loglikelihood of beta
    theta <- exp(beta[1]*z0 + beta[2]*z1 + beta[3]*z2)
    fbig <- f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    logfs <- log_fs(x,theta)
    d = fbig - alpha 
    
    lab = 0
    if(min(d)*max(d) <0 ){
      fl  <- function(lab){ sum(d/(1+lab*d))}  
      low= -1/max(d) + eeps
      up = -1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- sum(logfs) - sum(log( 1 + lab*d))
    -f
  }
  maxlike <- digam(n_big = n.big,n_small = n) + (n.big-n)*log(alpha)- optim(par = c(1,1,1),fn=prof.lik)$value 
  maxlike
}





##--------------------    conditional likelihood procedure     -----------------------#
##    case  3
max3.clike <-function(obs){ 
  n = nrow(obs)
  x = obs[, 1]
  y = obs[, 2]
  z0 = obs[,3]
  z1 = obs[,4]
  z2 = obs[,5]
  
  llc <- function(beta){
    theta <- exp(beta[1]*z0 + beta[2]*z1 + beta[3]*z2) 
    logfs <- log_fs(x,theta) 
    logfb <- log_fb(y,theta)
    -sum(logfs - logfb)
  }
  
  beta_tilde = optim(par = c(1,1,1), fn = llc)$par
  theta_tilde <- exp(beta_tilde[1]*z0+beta_tilde[2]*z1 + beta_tilde[3]*z2)
  
  N_tilde     = sum(1/(1- f_b(y, theta_tilde)))
  alpha_tilde = 1- n/N_tilde
  par_tilde   <- c(N_tilde, alpha_tilde, beta_tilde)
  return(par_tilde)
}


#####################################################################################################
###   case 3:   estimate the asymptotical variances of tilde(N), tilde(alpha) and tilde(beta0),tilde(beta1),tilde(beta2)

Variance3_est <- function(obs){
  x = obs[,1]
  y = obs[,2]
  z0 = obs[,3]
  z1 = obs[,4]
  z2 = obs[,5]
  qz = cbind(z0,z1,z2)
  n = nrow(obs)
  
  til_par = max3.clike(obs)
  til_N = til_par[1]
  til_a = til_par[2]
  til_b0 = til_par[3]
  til_b1 = til_par[4]
  til_b2 = til_par[5]
  
  til_theta = exp(til_b0*z0 + til_b1*z1 + til_b2*z2)
  
  f1.big <- f1_b(y,theta = til_theta)
  til_phi <- sum(exp( 2*til_theta*y))/til_N
  
  coef_V23 <- NULL 
  coef_V33 <- NULL #### calculate the coefficients of components of V_{33}^{-1}
  for (i in 1:n) {
    fn <- function(x){
      (1/til_theta[i] - 2*x +til_theta[i]*x^2)/exp(til_theta[i]*x)
    }
    int_val <- integrate(fn,lower = y[i],upper = Inf)[[1]]
    sec_val <- f1_b(y[i],til_theta[i])^2*exp(til_theta[i]*y[i])
    coef_V33[i] <- (int_val - sec_val )*til_theta[i]^2*exp(til_theta[i]*y[i])
    
    coef_V23 <- rbind(coef_V23,f1_b(y[i],til_theta[i])*til_theta[i]*exp(2*til_theta[i]*y[i])*qz[i,])
  }  
  til_V23 <-  -apply(coef_V23, 2, sum)/til_N
  V1_est  <-  -1/til_N*sum(coef_V33)
  Vz1_est  <- -1/til_N*sum(coef_V33*z1)
  Vz2_est  <- -1/til_N*sum(coef_V33*z2)## V_{z2}
  Vz1s_est  <- -1/til_N*sum(coef_V33*z1^2) ## V_{z1^2}
  Vz2s_est  <- -1/til_N*sum(coef_V33*z2^2) 
  Vz12_est  <- -1/til_N*sum(coef_V33*z1*z2)
  
  
  ############ the  components of V_{33}^{-1}  ########################
  Vs_est = Vz1_est*Vz2s_est - Vz2_est*Vz12_est
  Vt_est = Vz1s_est*Vz2s_est - Vz12_est^2
  Vk_est = V1_est*Vz2s_est - Vz2_est^2
  Vl_est = Vk_est*Vt_est - Vs_est^2
  
  inV_11 <-   Vz2s_est*Vt_est/Vl_est
  inV_12 <-  -Vz2s_est*Vs_est/Vl_est
  inV_13 <-  (Vs_est*Vz12_est - Vt_est*Vz2_est  )/Vl_est
  inV_22 <-   Vz2s_est*Vk_est/Vl_est
  inV_23 <-   (Vs_est*Vz2_est - Vk_est*Vz12_est  )/Vl_est
  inV_33 <-   (1 + (Vt_est*Vz2_est^2 - 2*Vs_est*Vz2_est*Vz12_est + Vk_est*Vz12_est^2 )/Vl_est   )/Vz2s_est
  
  com <- c(inV_11,inV_12,inV_13,inV_12,inV_22,inV_23,inV_13,inV_23,inV_33)
  til_inV33 <- matrix(com,ncol = 3,byrow = T,dimnames = NULL)
  
  Var_N <- til_phi - 1 - sum((til_V23 %*% til_inV33) *  til_V23)
  Var_A <- til_a*(1-til_a)
  Var_b0 <- -inV_11
  Var_b1 <- -inV_22
  Var_b2 <- -inV_33
  
  cbind(Var_N,Var_A,Var_b0,Var_b1,Var_b2)
}


