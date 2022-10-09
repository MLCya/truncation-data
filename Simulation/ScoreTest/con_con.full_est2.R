##  to calculate the type 1 error  and power  for qz=(1,z), z~N(1,1)     ###

library(multicool)
library(rootSolve)
library(cubature)
eeps <-  1.490116e-08





##----------      score test       ------------------------------------##

digam <- function(n_big, n_small){ 
  out=0
  if(n_small>0 &  n_big>=n_small)  out=sum(log(n_big-n_small+c(1:n_small)))    
  out
}




################################################################
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




##--------------------------------------------------------------------------## 
##  generate the quantiles of  "mixed" chi-squared distributions 
##    (mixtures of chi-square(n) and chi-square(n-1) 
##
##--------------------------------------------------------------------------##


qchibarsq <- function(p,df=1,mix=0.5) 
{
  n <- max(length(p),length(df),length(mix))
  df <- rep(df,length.out=n)
  mix <- rep(mix,length.out=n)
  p <- rep(p,length.out=n)
  tmpf2 <- function(p,df,mix) {
    if (df>1) {
      tmpf <- function(x) {
        pchibarsq(x,df,mix)-p
      }
      uniroot(tmpf,lower=qchisq(p,df-1),upper=qchisq(p,df))$root
    } else {
      newq <- (p-mix)/(1-mix)
      ifelse(newq<0,0,qchisq(newq,df=1))
    }
  }
  mapply(tmpf2,p,df,mix)
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






#--------------------------      generate data  ------------------------#
gen_obs <- function(beta,rate.G,tau,n_big){
  z0 <- rep(1,n_big)
  z <- rnorm(n_big,1,1)
  y <- rexp(n_big,rate = rate.G)
  nu <- rnorm(n_big,mean = 0, sd = sqrt(2))
  
  x <- NULL
  theta <- NULL
  for (i in 1:n_big) {
    hi <- exp(beta[1]*z0[i]+beta[2]*z[i] + nu[i]*sqrt(tau))
    x[i] <- rexp(1,rate = hi)
    theta[i] <- hi
  }
  
  ind <- (x>y)
  obs <- data.frame(obs.x = x[ind],obs.y = y[ind],z0.obs=z0[ind],z.obs=z[ind],theta.obs = theta[ind])
  return(obs)
}







### full likelihood esimation procedure     parametric model: theta = h(z,beta) = exp(beta0 +beta1*z)
maxlike.par <- function(obs){
  x <-  obs[,1]
  y <-  obs[,2]
  z0<-  obs[,3]
  z <-  obs[,4]
  n <- nrow(obs)
  
 
  prof.lik <- function(beta){  ###profile log likelihood of beta
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
  
 
  theta_Phat  <- exp(beta_Phat[1]*z0+beta_Phat[2]*z )  ##the eatimator of theta = h(z,beta)
  fbig_Phat   <- f_b(y,theta_Phat)
  fsmall_Phat <- f_s(x,theta_Phat)
  
  fbig.1_Phat <- min(fbig_Phat)
  fbig.n_Phat <- max(fbig_Phat)
  logfs_Phat <- log_fs(x,theta_Phat)
  
  fn.outer <- function(gam){ 
    alpha_ter <- 1/(1+gam^2) * fbig.1_Phat + gam^2/(1+gam^2)*fbig.n_Phat
    N     <-  MaxN.alpha(n_small = n,alpha = alpha_ter)
    d <- fbig_Phat  - alpha_ter
    
    lab = 0
    if(max(d)*min(d)<0){
      fl  <- function(lab){ sum(d/(1+lab*d))   }  
      low= -1/max(d) + eeps
      up = -1/min(d) - eeps 
      lab = uniroot(fl, c(low,up))$root
    }
    
    f <- digam(N,n) + (N-n)*log(alpha_ter) - sum(log(1+lab*d))##+sum(logfs_Phat)
    return(-f)
  }
  
  result <- optimize(fn.outer,interval = c(0,1))
  gama_Phat   <- result$minimum
  alpha_Phat <- 1/(1+gama_Phat^2)*fbig.1_Phat + gama_Phat^2/(1+gama_Phat^2)*fbig.n_Phat  ## the eatimator of alpha
  N_Phat      <- MaxN.alpha(n,alpha_Phat)       ## the eatimator of N
  
  d_Phat <- fbig_Phat  - alpha_Phat
  
  lab_Phat = 0
  if(max(d_Phat)*min(d_Phat)<0){
    fl  <- function(lab){ sum(d_Phat/(1+lab*d_Phat))   }  
    low= -1/max(d_Phat) + eeps
    up = -1/min(d_Phat) - eeps 
    lab_Phat = uniroot(fl, c(low,up))$root
  }
  
  par_Phat  <- c(N_Phat, alpha_Phat,beta_Phat,lab_Phat)
  return(par_Phat)
}

### conditional  likelihood esimation procedure
max.clike<-function(obs){ 
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  
  llc <- function(beta){
    theta   <- exp(beta[1]*z0 + beta[2]*z ) 
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


##--------------    to estimate c(beta)  --------------#
cb_est <- function(obs){
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  
  hat_par <- maxlike.par(obs)
  hat_N   <- hat_par[1]
  hat_alpha <- hat_par[2]
  hat_beta <- matrix(hat_par[3:4],ncol = 2,dimnames = NULL)
  hat_lab  <-  hat_par[5]
  
  hat_theta <- NULL
  hat_com <- NULL
  for (i in 1:n) {
    qz_i <- c(z0[i], z[i])
    hat_theta[i] <-  exp(hat_beta[1]*z0[i]+hat_beta[2]*z[i] )
    hat_pi <- 1/n * 1/(1+hat_lab*(f_b(y[i], hat_theta[i]) - hat_alpha )   )
    fn <- function(x){
      (f2_s(x,hat_theta[i])*hat_theta[i]+f1_s(x,hat_theta[i]))*(1/hat_theta[i] - x)
    }
    
    int_val <- integrate(fn,lower = y[i],upper = Inf)[[1]]
    sec_val <- (f2_b(y[i],hat_theta[i])*hat_theta[i]+f1_b(y[i],hat_theta[i]))*f1_b(y[i],hat_theta[i])*exp(hat_theta[i]*y[i])
    hat_com <- rbind(hat_com,  (int_val - sec_val)*hat_theta[i]^2*hat_pi*qz_i)
  }
  -1*apply(hat_com, 2, sum)
}

##----------------  to estimate V_{33}^{-1} when q(z) = (1,z1,z2)    ----------------##
inV33_est <- function(obs){
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  
  hat_par <- maxlike.par(obs)
  hat_N   <- hat_par[1]
  hat_alpha <- hat_par[2]
  hat_beta <- matrix(hat_par[3:4],ncol = 2,dimnames = NULL)
  hat_lab  <-  hat_par[5]
  
  hat_theta <- NULL  ## h(z,til_beta)
  hat_coef <- NULL
  for (i in 1:n) {
    hat_theta[i] <- exp(hat_beta[1]*z0[i]+hat_beta[2]*z[i])
    hat_pi <- 1/n * 1/(1+hat_lab*(f_b(y[i], hat_theta[i]) - hat_alpha )   )
    
    fn <- function(x){
      (1/hat_theta[i] - 2*x +hat_theta[i]*x^2)/exp(hat_theta[i]*x)
    }
    int_val <- integrate(fn,lower = y[i],upper = Inf)[[1]]
    sec_val <- f1_b(y[i],hat_theta[i])^2*exp(hat_theta[i]*y[i])
    hat_coef[i] <- (int_val - sec_val )*hat_theta[i]^2*hat_pi
  }  
  
  V1_est  <- -1*sum(hat_coef)
  Vz_est <-  -1*sum(hat_coef*z)
  Vzs_est <- -1*sum(hat_coef*z^2)## V_{z^2}
  
  V_d <- V1_est*Vzs_est - Vz_est^2
  
  ############ the  components of V_{33}^{-1}  ########################
  inV_11 <- Vzs_est/V_d
  inV_12 <- - Vz_est/V_d
  inV_22 <- V1_est/V_d

  com <- c(inV_11,inV_12,inV_12,inV_22)
  matrix(com,ncol = 2,byrow = T,dimnames = NULL)
}


###--------------- to estimate sigam_s^2 -----------------------####
sigmas2_est <- function(obs){
  n = nrow(obs)
  x  <- obs[,1]
  y  <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]

  hat_cb <- cb_est(obs)
  hat_inV33 <- inV33_est(obs)
  cv   <- c(hat_cb[1]*hat_inV33[1,1] + hat_cb[2]*hat_inV33[2,1],
            hat_cb[1]*hat_inV33[1,2] + hat_cb[2]*hat_inV33[2,2])
    
  hat_par <- maxlike.par(obs)
  hat_N   <-   hat_par[1]
  hat_alpha <- hat_par[2]
  hat_beta <- matrix(hat_par[3:4],ncol = 2,dimnames = NULL)
  hat_lab  <-  hat_par[5]
  
  
  hat_theta <- NULL  ## h(z,til_beta)
  hat_coef <- NULL
  for (i in 1:n) {
    #qz_i <- c(z0[i],z[i])
    hat_theta[i] <- exp(hat_beta[1]*z0[i] + hat_beta[2]*z[i])
    hat_pi <- 1/n * 1/(1+hat_lab*(f_b(y[i], hat_theta[i]) - hat_alpha )   )
    cvq <- cv[1] + cv[2]*z[i]
    #cvq <- sum( cv * qz_i)
    fn <- function(x){
      ((hat_theta[i]*x-2)*hat_theta[i]*x+(1-hat_theta[i]*x)*(1-cvq))^2/(hat_theta[i]*exp(hat_theta[i]*x))
    }
    
    int_val <- integrate(fn,lower = y[i],upper = Inf)[[1]]
    sec_val <- (f2_b(y[i],hat_theta[i])*hat_theta[i]+f1_b(y[i],hat_theta[i])*(1-cvq))^2 *exp(hat_theta[i]*y[i])  
    hat_coef[i] <- (int_val - sec_val )*hat_theta[i]^2 * hat_pi
  }
  sum(hat_coef)
}

###---------------     the standard score statistics     ---------------####
U_ts <- function(obs){
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z0 <- obs[,3]
  z <- obs[,4]
  
  hat_par <- maxlike.par(obs)
  hat_N   <- hat_par[1]
  hat_alpha <- hat_par[2]
  hat_beta <- matrix(hat_par[3:4],ncol = 2,dimnames = NULL)
  hat_lab  <-  hat_par[5]
  hat_sigms <- sqrt(sigmas2_est(obs)) 
  
  hat_theta <- NULL  ## h(z,til_beta)
  terp <- NULL
  for (i in 1:n) {
    hat_theta[i] <- exp(hat_beta[1]*z0[i]+ hat_beta[2]*z[i])
    fir_val <- hat_theta[i]*x[i]*(hat_theta[i]*x[i] -2) + (1- hat_theta[i]*x[i])
    sec_val <- -hat_theta[i]^2*y[i]^2 + hat_theta[i]*y[i]
    terp[i] <- fir_val + sec_val
  }
  1/sqrt(hat_N)*sum(terp)/hat_sigms
}


####     to simulate for  q(z) = (1,z)         ####
beta   <- c(-0.2,0.5)
rate.G <- 9
n_big  <- 500
nrep   <- 3000
K   <- 5
tau <- 0.02*c(0:K)  ## null&alternative
a   <- c(0.01,0.05)  ##significant  level

ty1val_uts <- NULL ##under H0, the value of the standard score statistic
pw1val_uts <- NULL 
pw2val_uts <- NULL 
pw3val_uts <- NULL 
pw4val_uts <- NULL 
pw5val_uts <- NULL 


ty1val_mchi <- NULL  ##under H0,the value of the mixed chi-square statistic
pw1val_mchi <- NULL
pw2val_mchi <- NULL
pw3val_mchi <- NULL
pw4val_mchi <- NULL
pw5val_mchi <- NULL


ty1_par_est <- NULL
pw1_par_est <- NULL
pw2_par_est <- NULL
pw3_par_est <- NULL
pw4_par_est <- NULL
pw5_par_est <- NULL

for (i in 1: nrep){
  ty1_obs  <- gen_obs(beta,rate.G,tau[1],n_big)
  pw1_obs  <- gen_obs(beta,rate.G,tau[2],n_big)
  pw2_obs  <- gen_obs(beta,rate.G,tau[3],n_big)
  pw3_obs  <- gen_obs(beta,rate.G,tau[4],n_big)
  pw4_obs  <- gen_obs(beta,rate.G,tau[5],n_big)
  pw5_obs  <- gen_obs(beta,rate.G,tau[6],n_big)
  
  ty1_uts <- U_ts(ty1_obs)
  pw1_uts <- U_ts(pw1_obs)
  pw2_uts <- U_ts(pw2_obs)
  pw3_uts <- U_ts(pw3_obs)
  pw4_uts <- U_ts(pw4_obs)
  pw5_uts <- U_ts(pw5_obs); 
  
  ty1val_uts <- c(ty1val_uts, ty1_uts)
  pw1val_uts <- c(pw1val_uts, pw1_uts)
  pw2val_uts <- c(pw2val_uts, pw2_uts)
  pw3val_uts <- c(pw3val_uts, pw3_uts)
  pw4val_uts <- c(pw4val_uts, pw4_uts)
  pw5val_uts <- c(pw5val_uts, pw5_uts)
  
  ty1val_mchi <- c(ty1val_mchi, max(ty1_uts,0)^2 )
  pw1val_mchi <- c(pw1val_mchi, max(pw1_uts,0)^2 )
  pw2val_mchi <- c(pw2val_mchi, max(pw2_uts,0)^2 )
  pw3val_mchi <- c(pw3val_mchi, max(pw3_uts,0)^2 )
  pw4val_mchi <- c(pw4val_mchi, max(pw4_uts,0)^2 )
  pw5val_mchi <- c(pw5val_mchi, max(pw5_uts,0)^2 )
  
  
  ty1_par_est <- rbind(ty1_par_est,c(maxlike.par(ty1_obs)[-6],max.clike(ty1_obs)) )
  pw1_par_est <- rbind(pw1_par_est,c(maxlike.par(pw1_obs)[-6],max.clike(pw1_obs)) )
  pw2_par_est <- rbind(pw2_par_est,c(maxlike.par(pw2_obs)[-6],max.clike(pw2_obs)) )
  pw3_par_est <- rbind(pw3_par_est,c(maxlike.par(pw3_obs)[-6],max.clike(pw3_obs)) )
  pw4_par_est <- rbind(pw4_par_est,c(maxlike.par(pw4_obs)[-6],max.clike(pw4_obs)) )
  pw5_par_est <- rbind(pw5_par_est,c(maxlike.par(pw5_obs)[-6],max.clike(pw5_obs)) )
  print(i)
}



ty1_indx <- cbind( (ty1val_mchi > qchibarsq(1-a[1],df=1)), (ty1val_mchi > qchibarsq(1-a[2],df=1))      )
pw1_indx <- cbind( (pw1val_mchi > qchibarsq(1-a[1],df=1)), (pw1val_mchi > qchibarsq(1-a[2],df=1))      )
pw2_indx <- cbind(  (pw2val_mchi > qchibarsq(1-a[1],df=1)), (pw2val_mchi > qchibarsq(1-a[2],df=1))      )
pw3_indx <- cbind( (pw3val_mchi > qchibarsq(1-a[1],df=1)), (pw3val_mchi > qchibarsq(1-a[2],df=1))      )
pw4_indx <- cbind( (pw4val_mchi > qchibarsq(1-a[1],df=1)), (pw4val_mchi > qchibarsq(1-a[2],df=1))      )
pw5_indx <- cbind( (pw5val_mchi > qchibarsq(1-a[1],df=1)), (pw5val_mchi > qchibarsq(1-a[2],df=1))      )


ty1_est <- apply(ty1_par_est , 2, mean)  
pw1_est <- apply(pw1_par_est , 2, mean)  
pw2_est <- apply(pw2_par_est , 2, mean)  
pw3_est <- apply(pw3_par_est , 2, mean)  
pw4_est <- apply(pw4_par_est , 2, mean)
pw5_est <- apply(pw5_par_est , 2, mean)


ty1_error <- apply(ty1_indx , 2, sum)/nrep
pw1 <- apply(pw1_indx, 2, sum)/nrep
pw2 <- apply(pw2_indx, 2, sum)/nrep
pw3 <- apply(pw3_indx, 2, sum)/nrep
pw4 <- apply(pw4_indx, 2, sum)/nrep
pw5 <- apply(pw5_indx, 2, sum)/nrep

val_uts  <- cbind(ty1val_uts, pw1val_uts,  pw2val_uts,  pw3val_uts,  pw4val_uts,  pw5val_uts)
val_mchi <- cbind(ty1val_mchi,pw1val_mchi, pw2val_mchi, pw3val_mchi, pw4val_mchi, pw5val_mchi  )
power_case2 <- data.frame(cbind(ty1_error,pw1,pw2,pw3,pw4,pw5), row.names = c("0.1%","1%","5%"))



