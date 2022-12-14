


############ supp-additional  scenario：F~weibull , G ~ exponentioal , z ~ U(1,10)
#####  qz = (1,z), theta = exp(beta[1] + beta[2]*z)
#####  case1：beta = (0, 0),  rate.G = 0.3, O_P1 =  0.7692307
#####  case2:  beta = (-0.5,0.3), rate.G = 0.1, O_P2 = 0.7340651
##### modelling: theta = h(z,beta) = exp(beta0 + beta1*z)


library(rootSolve)
library(cubature)

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

###### given alpha to estimate N
MaxN.alpha <- function(n_small, alpha){  
  myfun <- function(N){  
    tmp = N- n_small + (1:n_small)
    sum(1/tmp)+log_n(1 - alpha,n_small)  
  }
  multiroot(f = myfun,start = c(n_small),maxiter = 5000)$root
  #lower =  n_small
  #upper =  -n_small/log(alpha)  + n_small - 1
  #uniroot(myfun, c(lower+2, upper),extendInt = "yes",maxiter = 50000)$root ## ,extendInt = "yes"
}




f_s <- function(x,scale,shape = 1){
  shape/scale * (x/scale)^(shape -1)/exp((x/scale)^shape)
}


f_b <- function(x,scale,shape = 1){
  1 - exp(-(x/scale)^(shape))
}

f1_s <- function(x,scale,shape = 1){
  shape^2*x^(shape-1)*scale^(-shape-1)*((x/scale)^shape - 1)/exp((x/scale)^shape)
}

f2_s <- function(x,scale,shape=1){
  shape^2*x^(shape-1)/scale^(shape+2)/exp((x/scale)^shape)*(shape*(x/scale)^(2*shape) - 4*shape*(x/scale)^shape + shape +1)
  
}

f1_b <- function(x,scale, shape=1){
  -shape*x^shape*scale^(-shape-1)/exp((x/scale)^shape)
}

f2_b <- function(x,scale,shape=1){
  shape*x^shape/scale^(shape+2)/exp((x/scale)^shape)*(shape+1-shape*(x/scale)^shape)
}

log_fs <- function(x, scale, shape = 1){
  log(shape) - shape*log(scale) + (shape - 1)*log(x) - (x/scale)^shape
}

log_fb <- function(x, scale, shape = 1){
  log( 1 - exp(- (x/scale)^shape)  )
}




##---------  Calculate the sampling probability  ------------###

O_P <- function(beta,rate.G){# beta is vector of 2 dimension
  comf <- function(yz){#yz is vector
    y <- yz[1]
    z <- yz[2]
    scale <- exp(beta[1]+beta[2]*z)
    
    f_b(y,scale = scale)*rate.G*exp(-rate.G*y)*dunif(z,min = 1,max = 10)
  }
  adaptIntegrate(comf,lowerLimit = c(0,1),upperLimit = c(Inf,10),tol = 1e-05)$integral
}




##    generating  data: theta = h(z,beta) = exp(beta0 + beta1*z)
gen_obs <- function(beta,rate.G,n.big){#beta is a vector of 2 dimension,qz=(1,z)
  y = rexp(n.big,rate = rate.G)
  z = runif(n.big,1,10)
  scale_t = exp(beta[1] + beta[2]*z )
  x <- NULL
  for(i in 1:n.big){
    x[i] = rweibull(1,scale = scale_t[i],shape = 1)
  }
  
  ind = (x < y) ## right truncation
  obs <- data.frame(obs.x= x[ind],obs.y=y[ind],obs.z=z[ind],scale.obs = scale_t[ind])
  return(obs)
}









###  full likelihood procedure     parametric model:  theta = h(z,beta) = exp(beta0 + beta1*z)
maxlike.par <- function(obs){
  n <- nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z <- obs[,3]
  
  
  
  prof.lik <- function(beta){  ###profile loglikelihood of beta
    scale   <- exp(beta[1] + beta[2]*z) 
    fbig    <- f_b(y,scale = scale)
    fsmall  <- f_s(x,scale = scale)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,scale = scale) ## log(f(x,theta))
    
    
    fn.iner <- function(gam){ ## gam is a reparameter 
      alpha <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      N     <- MaxN.alpha(n_small = n,alpha = alpha )
      
      d <- fbig - alpha
      lab = 0  
      if(max(d) * min(d) < 0){
        fl  <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      
      fn <- digam(N,n) + (N-n)*log(1-alpha) - sum(log(1 + lab*d))
      return(-fn )
    }
    optimize(fn.iner,c(0,100))$objective - sum(logfs)
  }
  
  ## the eatimator of beta 
  out       <-  optim(c(0.5,0.5), fn=prof.lik)
  beta_Phat <-  out$par
  
  
  
  scale_Phat  <- exp(beta_Phat[1] + beta_Phat[2]*z )  ##the eatimator of theta = h(z,beta)
  fbig_Phat   <- f_b(y,scale = scale_Phat)
  fsmall_Phat <- f_s(x,scale = scale_Phat)
  
  fbig.1_Phat <- min(fbig_Phat)
  fbig.n_Phat <- max(fbig_Phat)
  logfs_Phat  <- log_fs(x,scale = scale_Phat)
  
  fn.outer <- function(gam){ ## ng = c(n.big,gama) is parametric vector
    alpha_ter <- 1/(1+gam^2) * fbig.1_Phat + gam^2/(1+gam^2)*fbig.n_Phat
    N <- MaxN.alpha(n_small = n, alpha = alpha_ter)
    d <- fbig_Phat  - alpha_ter
    
    lab = 0
    if(max(d)*min(d)<0){
      fl  <- function(lab){ sum(d/(1+lab*d))   }  
      low= -1/max(d) + eeps
      up = -1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- digam(N,n) + (N-n)*log(1-alpha_ter) +sum(logfs_Phat) - sum(log_n(1+lab*d,n))
    return(-f)
  }
  
  result <- optimize(fn.outer,interval = c(0,100))
  gama_Phat   <- result$minimum
  alpha_Phat <- 1/(1+gama_Phat^2)*fbig.1_Phat + gama_Phat^2/(1+gama_Phat^2)*fbig.n_Phat  ## the eatimator of alpha
  
  N_Phat  <- MaxN.alpha(n,alpha_Phat)       ## the eatimator of N
  #theta_Phat.m <- mean(theta_Phat)
  par_Phat  <- c(N_Phat, alpha_Phat,beta_Phat)
  maxlike   =  -out$value 
  list(par_Phat,maxlike)
}




##------------- calculating likelihood ratio for parameters  ------------------------------##
###   1、N=N0 is known 

maxlike.N <- function(obs,n.big){
  x <- obs[,1]
  y <- obs[,2]
  z <- obs[,3]
  n <- nrow(obs)
  N <- n.big
  
  prof.lik <- function(beta){#  profile loglikelihood of beta
    scale <- exp(beta[1] + beta[2]*z)
    fbig  <- f_b(y,scale = scale)
    fsmall <- f_s(x,scale = scale)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,scale = scale)
    
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
      fn <- (N - n)*log(1 - alpha) - sum(log_n(1+lab*d,n))
      -fn
    }
    optimize(fn.inner,interval = c(0,100))$objective - sum(logfs )
  }
  maxlike <- digam(N,n) - optim(par = c(0.5,0.5),fn=prof.lik)$value  
  maxlike
}



###   2、 alpha = alpha0  is known

maxlike.A  <- function(obs,alpha){
  x <-  obs[,1]
  y <-  obs[,2]
  z  <- obs[,3]
  n <- nrow(obs)
  alpha <- alpha
  n.big <-  MaxN.alpha(n, alpha = alpha) ## the eatimator of N
  
  prof.lik  <- function(beta){  ## profile loglikelihood of beta
    scale  <- exp(beta[1] + beta[2]*z)
    fbig   <-  f_b(y,scale = scale)
    fsmall <-  f_s(x,scale = scale)
    
    logfs <- log_fs(x,scale = scale)
    d = fbig - alpha 
    
    lab = 0
    if(min(d)*max(d) <0 ){
      fl  <- function(lab){ sum(d/(1+lab*d))}  
      low= -1/max(d) + eeps
      up = -1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f  <-  sum(logfs) - sum(log( 1 + lab*d))
    -f
  }
  maxlike <- digam(n_big = n.big,n_small = n) + (n.big - n)*log(1-alpha) - optim(par = c(0.1,0.1),fn=prof.lik)$value 
  maxlike
}


### 3、  beta = beta0  is known


maxlike.B0 <- function(obs,b0){##given  beta0
  x <- obs[,1]
  y <- obs[,2]
  z <- obs[,3]
  n <- nrow(obs)
  
  pB1 <- function(b1){
    beta  <-  c(b0,b1)
    scale <- exp(beta[1]+beta[2]*z)  # true theta = h(z)
    fbig   <-  f_b(y,scale=scale)
    fsmall <-  f_s(x, scale = scale)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x,scale = scale)
    
    fng <- function(gam){ ## ng = c(n.big,gama) is parametric vector
      alpha_ter <- 1/(1+gam^2) * fbig.1 +  gam^2/(1+gam^2)*fbig.n
      n.big <- MaxN.alpha(n_small = n, alpha = alpha_ter)
      
      d <- fbig - alpha_ter
      
      lab = 0
      if(min(d)*max(d) <0 ){
        fl <- function(lab){ sum(d/(1+lab*d))}  
        low= -1/max(d) + eeps
        up = -1/min(d) - eeps 
        lab = uniroot(fl, c(low,up))$root
      }
      f <- digam(n.big, n) + (n.big - n)*log(1-alpha_ter) + sum(logfs) - sum(log(1 + lab*d)) 
      -f 
    }
    optimize(fng,c(0,100))$objective
  }
  
  maxlike  <- -1*optimize(pB1,c(-10,10))$objective
  maxlike  
}


maxlike.B1 <- function(obs,b1){##given  beta1
  x <- obs[,1]
  y <- obs[,2]
  z <- obs[,3]
  n <- nrow(obs)
  
  pB0 <- function(b0){
    beta  <-  c(b0,b1)
    scale <-  exp(beta[1]+beta[2]*z)  # true theta = h(z)
    fbig   <-  f_b(y,scale = scale)
    fsmall <-  f_s(x, scale = scale)
    
    fbig.1 <- min(fbig)
    fbig.n <- max(fbig)
    logfs <- log_fs(x, scale = scale)
    
    fn <- function(gam){ ## ng = c(n.big,gama) is parametric vector
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
      f <- digam(n.big, n) + (n.big - n)*log(1- alpha_ter) + sum(logfs) - sum(log(1 + lab*d)) 
      -f 
    }
    optimize(fn,c(0,100))$objective
  }
  
  maxlike  <- -1*optimize(pB0,c(-10,10))$objective
  maxlike  
}



###      conditional likelihood  estimating procedure: theta = h(z,beta) = exp(beta0 + beta1*z)   ###########
max.clike<-function(obs){ 
  n = nrow(obs)
  x = obs[,1]
  y = obs[,2]
  z = obs[,3]
  
  llc <- function(beta){
    scale <- exp(beta[1] + beta[2]*z ) 
    logfs <- log_fs(x,scale = scale) 
    logfb <- log_fb(y,scale = scale)
    -sum(logfs - logfb)
  }
  
  beta_tilde = optim(par = c(0.5,0.5),fn = llc)$par
  scale_tilde <- exp(beta_tilde[1]+beta_tilde[2]*z)
  
  N_tilde     = sum(1/(f_b(y, scale_tilde)))
  alpha_tilde = n/N_tilde
  par_tilde   <- c(N_tilde, alpha_tilde, beta_tilde)
  return(par_tilde)
}




#####################################################################################################
##supp-additional case:    estimate the asymptotical variances of tilde(N), tilde(alpha) and tilde(beta0),tilde(beta1).
Variance_est <- function(obs,shape = 1){
  x = obs[,1]
  y = obs[,2]
  z = obs[,3]
  
  qz = cbind(1,z)
  n = nrow(obs)
  
  til_par = max.clike(obs)
  til_N = til_par[1]
  til_a = til_par[2]
  til_b0 = til_par[3]
  til_b1 = til_par[4]
  
  
  til_scale = exp(til_b0 + til_b1*z)
  
  f1.big <- f1_b(y,scale = til_scale)
  f.big <-  f_b(y,scale = til_scale)
  til_phi <- sum(1/f.big^2)/til_N
  
  coef_V23 <- NULL 
  coef_V33 <- NULL ## ## calculate the coefficients of components of V_{33}^{-1}
  for (i in 1:n) {
    fn <- function(x){
      #shape*til_scale[i]*((x/til_scale[i])^shape - 1)
      f1_s(x,scale = til_scale[i])^2*til_scale[i]^2/f_s(x,scale = til_scale[i])
    }
    
    int_val <- integrate(fn,lower = 0,upper = y[i])[[1]]
    sec_val <- f1.big[i]^2*til_scale[i]^2/f.big[i]
    #sec_val <- f1_b(y[i], scale = til_scale[i])^2*til_scale[i]^2*exp( (x/til_scale[i])^shape )/(exp( (x/til_scale[i])^shape ) - 1 )
    
    coef_V33  <- rbind(coef_V33, (int_val - sec_val )/f.big[i]/til_N)
    coef_V23  <- rbind(coef_V23, f1_b(y[i],scale = til_scale[i])*til_scale[i]*qz[i,]/f.big[i]^2)
  }
  
  til_V23 <-  -apply(coef_V23, 2, sum)/til_N
  V1_est  <- -1*sum(coef_V33)
  Vz_est  <-  -1*sum(coef_V33*z)
  Vzs_est <- -1*sum(coef_V33*z^2)## V_{z^2}
  V_d <- V1_est*Vzs_est - Vz_est^2
  
  ############ the   components of V_{33}^{-1}  ########################
  inV_11 <- Vzs_est/V_d
  inV_12 <- -1*Vz_est/V_d
  inV_22 <- V1_est/V_d
  com <- c(inV_11,inV_12,inV_12,inV_22)
  
  til_inV33 <- matrix(com,ncol = 2,byrow = T,dimnames = NULL)
  
  Var_N <- til_phi - 1 - sum((til_V23 %*% til_inV33) *  til_V23)
  
  Var_A <- til_a*(1-til_a)
  Var_b0 <- -1*inV_11
  Var_b1 <- -1*inV_22
  
  cbind(Var_N,Var_A,Var_b0,Var_b1)
}

