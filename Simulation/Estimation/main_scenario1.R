
############    main scenrio 1：F&G exponential distribution,  z ~ N(1,2)
####### qz = (1,z), theta = exp(beta[1] + beta[2]*z)
#####  case1： beta = (0.2, 0),  rate.G = 5, O_P1 = 0.8036772
#####  case2:  beta = (0.3,0.5), rate.G = 9, O_P2 = 0.7802803





library(rootSolve)
library(cubature)
library(latex2exp)
library(kedd)
eeps <-  1.490116e-08 

digam <- function(n_big, n_small){ 
  out=0
  if(n_small>0 &  n_big>=n_small)  out=sum(log(n_big-n_small+c(1:n_small)))    
  out
}

#####  define a new log function  log_n ###

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

######given alpha to estimate N
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


## both F&G  exponential distribution
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
O_P <- function(beta,rate.G){# beta is vector of 2 dimension
  comf <- function(yz){#yz is vector
    theta <- exp(beta[1]+beta[2]*yz[2])
    #pexp(yz[1],theta)*dexp(yz[1],rate.G)*dnorm(yz[2],0,1)
    f_b(yz[1],theta = theta)*f_s(yz[1],rate.G)*dnorm(yz[2],mean = 1,sd = sqrt(2))
  }
  1 - adaptIntegrate(comf,lowerLimit = c(0,-Inf),upperLimit = c(Inf,Inf),tol = 1e-05)$integral
}



##          generating  data: theta = h(z,beta) = exp(beta0 + beta1*z)
gen_obs <- function(beta,rate.G,n.big){#beta is a vector of 2 dimension,qz=(1,z)
  y = rexp(n.big,rate = rate.G)
  z0 = rep(1,n.big)
  z = rnorm(n.big,1,1)
  theta_t = exp(beta[1]+beta[2]*z)
  x = rexp(n.big,rate = theta_t)
  
  ind = (x > y)
  obs <- data.frame(obs.x= x[ind],obs.y=y[ind],obs.z0=z0[ind],obs.z=z[ind],theta.obs = theta_t[ind])
  return(obs)
}




###  full likelihood procedure     parametric model: theta = h(z,beta) = exp(beta0 + beta1*z)
maxlike.par <- function(obs){
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
  
  par_Phat  <- c(N_Phat, alpha_Phat,beta_Phat)
  maxlike   =  -out$value 
  list(par_Phat,maxlike)
}

##------------- calculating likelihood ratio for parameters  ------------------------------##
###   1、N=N0 is known 
maxlike.N <- function(obs,n.big){
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



###   2、 alpha = alpha0  is known

maxlike.A  <- function(obs,alpha){
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




### 3、  beta = beta0  is known


maxlike.B0 <- function(obs,b0){##given beta0
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z  <- obs[,4]
  n <- nrow(obs)
  
  pB1 <- function(b1){
    beta  <-  c(b0,b1)
    theta <- exp(beta[1]*z0+beta[2]*z)  # true theta = h(z)
    fbig   <-  f_b(y,theta)
    fsmall <- f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta)
    
    fng <- function(gam){ ## ng = c(n.big,gama) is parametric vector
      ## reparametrics
      alpha_ter <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      n.big <- MaxN.alpha(n_small = n,alpha = alpha_ter )
      
      d <- fbig - alpha_ter
      
      lab = 0
      if(min(d)*max(d) <0 ){
        fl <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      f <- digam(n.big, n) + (n.big - n)*log(alpha_ter) - sum(log(1 + lab*d)) 
      -f 
    }
    -sum(logfs) + optimize(fng,c(0,1))$objective
  }
  
  maxlike  <- -optimize(pB1,c(-10,10))$objective
  maxlike  
}


maxlike.B1 <- function(obs,b1){##given beta1
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z  <- obs[,4]
  n <- nrow(obs)
  
  pB0 <- function(b0){
    beta  <-  c(b0,b1)
    theta <-  exp(beta[1]*z0+beta[2]*z)  # true theta = h(z)
    fbig   <-  f_b(y,theta)
    fsmall <-  f_s(x,theta)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,theta)
    
    fn <- function(gam){ ## ng = c(n.big,gama) is parametric vector
      ## reparametrics
      alpha_ter <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      n.big <- MaxN.alpha(n_small = n,alpha = alpha_ter )
      
      d <- fbig - alpha_ter
      
      lab = 0
      if(min(d)*max(d) <0 ){
        fl <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + 10*eeps
        up = -1/min(d) - 10*eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      f <- digam(n.big, n) + (n.big - n)*log(alpha_ter) - sum(log(1 + lab*d)) 
      -f 
    }
    -sum(logfs) + optimize(fn,c(0,1))$objective
  }
  
  maxlike  <-   -optimize(pB0,c(-10,10))$objective
  maxlike  
}




###      conditional likelihood  estimating procedure : theta = h(z,beta) = exp(beta0 + beta1*z)
max.clike<-function(obs){
  n = nrow(obs)
  x = obs[, 1]
  y = obs[, 2]
  z0 = obs[, 3]
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



### estimate the asymptotical variances of tilde(N), tilde(alpha) and tilde(beta0),tilde(beta1).
Variance_est <- function(obs){
  x = obs[,1]
  y = obs[,2]
  z0 = obs[,3]
  z = obs[,4]
  qz = cbind(z0,z)
  n = nrow(obs)
  
  til_par = max.clike(obs)
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
  
  ############ the   components of V_{33}^{-1}  ########################
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


