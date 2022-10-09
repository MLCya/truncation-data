
##-------------------  Final version of  Real Analysis : Rotterdam breast cancer data   ----------------------------#

#---------------------------------------------------------------------------#    
#                                                                       #
#                 Some useful function                                  #
#                                                                       #
#---------------------------------------------------------------------------#
library(survival)
library(rootSolve) 
eeps <-  1.490116e-08  ####precision


digam <- function(n_big, n_small){ 
  out=0
  if(n_small > 0 &  n_big>=n_small)  out=sum(log(n_big-n_small+c(1:n_small)))    
  out
}

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

###     estimate alpha givwn N estimator
MaxN.alpha <- function(n_small, alpha){  
  myfun <- function(N){  
    tmp = N +1 - (1:n_small)
    sum(1/tmp) +log_n(1- alpha,n_small)  
  }
  multiroot(f = myfun,start = c(n_small),maxiter = 5000)$root
  #lower =  n_small
  #upper =  -n_small/log(alpha)  + n_small - 1
  #uniroot(myfun, c(lower+2, upper),extendInt = "yes",maxiter = 50000)$root ## ,extendInt = "yes"
}



 
#### some usrful density or distribution function
f_b <- function(x,theta){
  1 - exp(-theta*x)
}


f_s <- function(x,theta){
  theta*exp(-theta*x) 
}

log_fs <- function(x,theta){
  log(theta) - theta*x
}

log_fb <- function(x,theta){## log(F(y,theta))
  log(1 - exp(-theta*x))
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



#----------------------------------------------------------------------------------------# 

#       parametric model:  theta = h(z,beta) = exp(beta^T q(z)), q(z) = (1,z1,z2,z3,z4,z5,z6,z7)

#       Full likelihood to estimate parameters
#---------------------------------------------------------------------------------------------#

maxlike.par <- function(obs,inti = rep(0,8)){
  n <- nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z1 <- obs[,3]
  z2 <- obs[,4]
  z3 <- obs[,5]
  z4 <- obs[,6]
  z5 <- obs[,7]
  z6 <- obs[,8]
  z7 <- obs[,9]
  
  

  
  prof.lik <- function(beta){  ###profile loglikelihood of beta
    theta   <- exp(beta[1] + beta[2]*z1 + beta[3]*z2 +  beta[4]*z3+  beta[5]*z4+  beta[6]*z5 +  beta[7]*z6 +  beta[8]*z7)
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
      
      fn <- digam(N,n) + (N-n)*log(1- alpha) - sum(log_n(1+lab*d,n))
      return(- fn )
    }
    optimize(fn.iner,c(0,1))$objective - sum(logfs)
  }
  
  
  
  ## the eatimator of beta 
  out <- optim(par = inti, fn=prof.lik)
  beta_Phat <-  out$par
  
  maxlike   <-  -out$value
  ############################################################################################# 
  
  theta_Phat  <- exp(beta_Phat[1] + beta_Phat[2]*z1 +beta_Phat[3]*z2+beta_Phat[4]*z3+beta_Phat[5]*z4+beta_Phat[6]*z5+beta_Phat[7]*z6+beta_Phat[8]*z7)  ##the eatimator of theta = h(z,beta)
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
    
    f <- digam(N,n) + (N-n)*log(1 - alpha_ter) - sum(log_n(1+lab*d,n))
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
  list(par_Phat,maxlike)
}

#----------------------------------------------------------------------------------------# 

#       parametric model:  theta = h(z,beta) = exp(beta^T q(z)), q(z) = (1,z1,z2,z3,z4,z5,z6,z7)

#       Condition  likelihood to estimate parameters
#---------------------------------------------------------------------------------------------#


max.clike<-function(obs,inti=rep(0,8)){ 
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z1 <- obs[,3]
  z2 <- obs[,4]
  z3 <- obs[,5]
  z4 <- obs[,6]
  z5 <- obs[,7]
  z6 <- obs[,8]
  z7 <- obs[,9]
  
  
  llc <- function(beta){
    theta   <- exp(beta[1] + beta[2]*z1 + beta[3]*z2 +  beta[4]*z3+  beta[5]*z4 + beta[6]*z5+ beta[7]*z6+ beta[8]*z7)
    logfs <- log_fs(x,theta) 
    logfb <- log_fb(y,theta)
    -sum(logfs - logfb)
  }
  out <- optim(par = inti, fn = llc)
  
  beta_tilde = out$par
  maxclike <- -1*out$value
  
  
  theta_tilde <- exp(beta_tilde[1] + beta_tilde[2]*z1 + beta_tilde[3]*z2 + beta_tilde[4]*z3 + beta_tilde[5]*z4+ beta_tilde[6]*z5+ beta_tilde[7]*z6+ beta_tilde[8]*z7)
  N_tilde     = sum(1/( f_b(y, theta_tilde)))
  alpha_tilde = n/N_tilde
  par_tilde   <- c(N_tilde, alpha_tilde, beta_tilde)
  list(par_tilde, maxclike)
}


#----------------------------------------------------------------------------------------#
#      asymptotic distribution & asymptotic covariance & (grouped) distribution and its confidence band
#      hat_par is Full likelihood estimates ,Siglev is confidence level
#----------------------------------------------------------------------------------------#




Asyvar.full <- function(obs,hat_par,Siglev){
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z1 <- obs[,3]
  z2 <- obs[,4]
  z3 <- obs[,5]##size1 = "<20"
  z4 <- obs[,6]##size2 = "20-50"
  z5 <- obs[,7]
  z6 <- obs[,8]
  z7 <- obs[,9]
  
  
  hat_N   <- hat_par[1]
  hat_alpha <- hat_par[2]
  hat_beta <- matrix(hat_par[3:10],ncol = 8,dimnames = NULL)
  hat_lab  <-  hat_par[11]
  
  hat_theta <- NULL
  hat_pi <- NULL
  hat_fbig <- NULL
  hat_f1.big <- NULL
  hat_firint <- NULL ## intergrated part
  
  
  V22_comp   <- NULL
  V23_comp <- NULL ## the components of V23
  V24_comp <- NULL
  V33_comp.coef <- NULL ## the coefficient of the components of V33
  V44_comp <- NULL
  
  for (i in 1:n) {
    hat_theta[i] <- exp(hat_beta[1] + hat_beta[2]*z1[i]+ hat_beta[3]*z2[i]+ hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                        + hat_beta[6]*z5[i]+ hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    hat_pi[i] <- 1/(1  + hat_lab*(f_b(y[i],hat_theta[i]) - hat_alpha))/n
    
    hat_fbig[i]  <- f_b(y[i],hat_theta[i])
    hat_f1.big[i] <- f1_b(y[i],hat_theta[i])
    
    
    fir_int <- function(t){
      f1_s(t,hat_theta[i])^2/f_s(t,hat_theta[i])
    }
    
    hat_firint[i] <- integrate(fir_int,lower = 0,upper = y[i])$value
    
    
    qz_i <- c(1,z1[i],z2[i],z3[i],z4[i],z5[i],z6[i],z7[i])
    
    V22_comp[i] <-  hat_pi[i]/hat_fbig[i] 
    
    V23_comp <- rbind(V23_comp, hat_f1.big[i]*hat_theta[i]*hat_pi[i]*qz_i/hat_fbig[i])
    
    V24_comp[i] <- hat_pi[i]/hat_fbig[i] 
    
    V33_comp.coef[i] <- ( hat_firint[i] - hat_f1.big[i]^2/hat_fbig[i]  )*hat_pi[i]*hat_theta[i]^2
    V44_comp[i] <- (hat_fbig[i] - hat_alpha)^2*hat_pi[i]/hat_fbig[i]
  }
  
  hat_V22 <- -1/(1-hat_alpha) + sum(V22_comp)
  hat_V23 <- -1*apply(V23_comp,2,sum)
  hat_V24 <- hat_alpha*sum(V24_comp)
  hat_V34 <- hat_alpha^2*hat_V23
  hat_V44 <- hat_alpha^2*sum(V44_comp)
  
  V33_11  <-  -1*sum(V33_comp.coef*1*1)
  V33_z1  <-  -1*sum(V33_comp.coef*z1) ##V_z1
  V33_z2  <-  -1*sum(V33_comp.coef*z2);V33_z3  <-  -1*sum(V33_comp.coef*z3);V33_z4  <-  -1*sum(V33_comp.coef*z4)
  V33_z5  <-  -1*sum(V33_comp.coef*z5);V33_z6  <-  -1*sum(V33_comp.coef*z6);V33_z7  <-  -1*sum(V33_comp.coef*z7)
  
  V33_z11  <-  -1*sum(V33_comp.coef*z1*z1) ##V_z11 = V_z1z1
  V33_z12  <-  -1*sum(V33_comp.coef*z1*z2);V33_z13  <-  -1*sum(V33_comp.coef*z1*z3);V33_z14  <-  -1*sum(V33_comp.coef*z1*z4)
  V33_z15  <-  -1*sum(V33_comp.coef*z1*z5);V33_z16  <-  -1*sum(V33_comp.coef*z1*z6);V33_z17  <-  -1*sum(V33_comp.coef*z1*z7)
  
  V33_z22  <-  -1*sum(V33_comp.coef*z2*z2) ##V_z22 = V_z2z2
  V33_z23  <-  -1*sum(V33_comp.coef*z2*z3);V33_z24  <-  -1*sum(V33_comp.coef*z2*z4)
  V33_z25  <-  -1*sum(V33_comp.coef*z2*z5);V33_z26  <-  -1*sum(V33_comp.coef*z2*z6);V33_z27  <-  -1*sum(V33_comp.coef*z2*z7)
  
  V33_z33  <-  -1*sum(V33_comp.coef*z3*z3);V33_z34  <-  -1*sum(V33_comp.coef*z3*z4)
  V33_z35  <-  -1*sum(V33_comp.coef*z3*z5);V33_z36  <-  -1*sum(V33_comp.coef*z3*z6);V33_z37  <-  -1*sum(V33_comp.coef*z3*z7)
  
  
  V33_z44  <-  -1*sum(V33_comp.coef*z4*z4)
  V33_z45  <-  -1*sum(V33_comp.coef*z4*z5);V33_z46  <-  -1*sum(V33_comp.coef*z4*z6);V33_z47  <-  -1*sum(V33_comp.coef*z4*z7)
  
  V33_z55  <-  -1*sum(V33_comp.coef*z5*z5);V33_z56  <-  -1*sum(V33_comp.coef*z5*z6);V33_z57  <-  -1*sum(V33_comp.coef*z5*z7)
  
  V33_z66  <-  -1*sum(V33_comp.coef*z6*z6);V33_z67  <-  -1*sum(V33_comp.coef*z6*z7)
  
  V33_z77  <-  -1*sum(V33_comp.coef*z7*z7)
  
  hat_V33 <- matrix(NA,8,8)
  top_tri <- matrix(c(V33_11,V33_z1,V33_z2,V33_z3,V33_z4,V33_z5,V33_z6,V33_z7,
                      NA,V33_z11,V33_z12,V33_z13,V33_z14,V33_z15,V33_z16,V33_z17,
                      NA,NA,V33_z22,V33_z23,V33_z24,V33_z25,V33_z26,V33_z27,
                      NA,NA,NA,V33_z33,V33_z34,V33_z35,V33_z36,V33_z37,
                      NA,NA,NA,NA,V33_z44, V33_z45, V33_z46, V33_z47,
                      NA,NA,NA,NA,NA,V33_z55,V33_z56,V33_z57,
                      NA,NA,NA,NA,NA,NA,V33_z66,V33_z67,
                      NA,NA,NA,NA,NA,NA,NA, V33_z77),8,8,byrow = TRUE)
  for(i in 1:8){
    for(j in 1:8){
      if(i>=j){
        hat_V33[i,j]<- top_tri[j,i]
      }else{
        hat_V33[i,j]<- top_tri[i,j]
      }
    }
  }
  
  hat_V33.invers <- solve(hat_V33)
  hat_phi <- sum(hat_pi/hat_fbig)
  
  N_asyVar <-c(  hat_phi - 1 - t(as.matrix(hat_V23)) %*%  hat_V33.invers %*% as.matrix(hat_V23)  ) 
  alpha_asyVar <- hat_alpha*(1-hat_alpha)
  beta_asyVar <-  -1* hat_V33.invers
  
  
  ## -----------    overall&subgrouped pointwise confidence band   -----------------------------------#
  
  Z_Sig <- qnorm(1-Siglev/2)  
  
  W22 <- -hat_V22 + hat_V24^2/hat_V44
  W23 <- -hat_V23 + hat_V34*hat_V24/hat_V44
  W33 <- -hat_V33 + as.matrix(hat_V34)%*% t(as.matrix(hat_V34))/hat_V44
  
  W_prime <- as.matrix( rbind(c(W22,W23),cbind(c(W23),W33)  ))
  W_prime.invers <- solve(W_prime)
  
  
  Delta1 <- function(t){
    fbig_t <- NULL
    fbig1_t <- NULL
    sec_coef <- NULL
    for (i in 1:n) {
      qz_i <- c(1,z1[i],z2[i],z3[i],z4[i],z5[i],z6[i],z7[i])
      
      fbig_t[i] <- f_b(t,hat_theta[i])
      fbig1_t[i] <- f1_b(t,hat_theta[i])
      
      sec_coef <- rbind(sec_coef, (fbig1_t[i] - fbig_t[i]* hat_f1.big[i]/hat_fbig[i]  )*hat_theta[i]*hat_pi[i]*qz_i )
      
    }
    
    fir_comp <- sum(fbig_t*hat_pi)/hat_alpha
    sec_comp <- apply(sec_coef, 2, sum)
    c(fir_comp, sec_comp)
  }

 #------------------------- overall&subgrouped distribution function -------------------------------------#
  idx_size1 <- (z3 ==1 &  z4 ==0)
  idx_size2 <- (z3 ==0 &  z4 ==1)
  idx_size3 <- (z3 ==0 &  z4 ==0)
  
  fz1_f.hat <- hat_pi*idx_size1/sum( hat_pi*idx_size1 ) ## size = "<20"density 
  fz2_f.hat <- hat_pi*idx_size2/sum( hat_pi*idx_size2 )
  fz3_f.hat <- hat_pi*idx_size3/sum( hat_pi*idx_size3 )
  
  
  fb.band <- function(t){
    fb_val <- NULL
    for (k in 1:n) {
      fb_val[k] <- f_b(t,hat_theta[k])
    }
    f_big <- sum(hat_pi   * fb_val )
    fb_val.size1 <- sum(fz1_f.hat   * fb_val )
    fb_val.size2 <- sum(fz2_f.hat   * fb_val )
    fb_val.size3 <- sum(fz3_f.hat   * fb_val )
    
    critcal <- sqrt(  t(as.matrix(Delta1(t)))%*% W_prime.invers %*% as.matrix(Delta1(t)) )/sqrt(hat_N)
    
    c(f_big, fb_val.size1, fb_val.size2,fb_val.size3, critcal)
    
  }

  rfst <- seq(0,20,by = 0.1)
  ter_val    <- NULL
  for (i in 1: length(rfst)) {
    ter_val <- rbind(ter_val,  fb.band(rfst[i]) )
  }
  
  fbig_int.or <- data.frame(fbig = ter_val[,1], fbig.u = ter_val[,1] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,1] -  Z_Sig* ter_val[,5])  ##overall
  
  fbig_int.s1 <- data.frame(fbig = ter_val[,2], fbig.u = ter_val[,2] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,2] -  Z_Sig* ter_val[,5])  ## Size1
  
  fbig_int.s2 <- data.frame(fbig = ter_val[,3], fbig.u = ter_val[,3] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,3] -  Z_Sig* ter_val[,5])  ## Size2
  fbig_int.s3 <- data.frame(fbig = ter_val[,4], fbig.u = ter_val[,4] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,4] -  Z_Sig* ter_val[,5])  ## Size3
  
  
  ##------------------- output  -------------------------#
  list("N_asyVar"= N_asyVar, "alpha_asyVar" = alpha_asyVar, "beta_asyVar"=beta_asyVar,"fbig_intval" =  fbig_int.or,
       "fbig_intval.size1" = fbig_int.s1,"fbig_intval.size2" = fbig_int.s2,"fbig_intval.size3" = fbig_int.s3  )
}





#----------------------------------------------------------------------------------------#

#              asymptotic distribution & asymptotic covariance & (grouped) distribution and its confidence band            

#              til_par is conditional likelihood estimates ,Siglev is confidence level

#----------------------------------------------------------------------------------------#

Asyvar.con <- function(obs,til_par,Siglev){
  n = nrow(obs)
  x <- obs[,1]
  y <- obs[,2]
  z1 <- obs[,3]
  z2 <- obs[,4]
  z3 <- obs[,5]##size1 = "<20"
  z4 <- obs[,6]##size2 = "20-50"
  z5 <- obs[,7]
  z6 <- obs[,8]
  z7 <- obs[,9]
  
  
  til_N   <- til_par[1]
  til_alpha <- til_par[2]
  til_beta <- matrix(til_par[3:10],ncol = 8,dimnames = NULL)
  
  
  til_theta <- NULL
  til_pi <- NULL
  til_fbig <- NULL
  til_f1.big <- NULL
  til_firint <- NULL 
  
  V23_comp <- NULL
  V33_comp.coef <- NULL 
  
  for (i in 1:n) {
    til_theta[i] <- exp(til_beta[1] + til_beta[2]*z1[i]+ til_beta[3]*z2[i]+ til_beta[4]*z3[i]+ til_beta[5]*z4[i]
                        + til_beta[6]*z5[i]+ til_beta[7]*z6[i]+ til_beta[8]*z7[i])
    
    til_fbig[i]  <- f_b(y[i],til_theta[i])
    til_f1.big[i] <- f1_b(y[i],til_theta[i])
    
    til_pi[i] <- 1/til_fbig[i]/til_N
    
    fir_int <- function(t){
      f1_s(t,til_theta[i])^2/f_s(t,til_theta[i])
    }
    
    til_firint[i] <- integrate(fir_int,lower = 0,upper = y[i])$value
    
    
    qz_i <- c(1,z1[i],z2[i],z3[i],z4[i],z5[i],z6[i],z7[i])
    
    V23_comp <- rbind(V23_comp, til_f1.big[i]* til_theta[i]*qz_i*til_pi[i]/til_fbig[i])
    
    V33_comp.coef[i] <- ( til_firint[i] - til_f1.big[i]^2/til_fbig[i]  )*til_pi[i]*til_theta[i]^2
  }
  
  til_V23 <- -1*apply(V23_comp,2,sum)
  
  V33_11  <-  -1*sum(V33_comp.coef*1*1)
  
  V33_z1  <-  -1*sum(V33_comp.coef*z1) ##V_z1
  V33_z2  <-  -1*sum(V33_comp.coef*z2);V33_z3  <-  -1*sum(V33_comp.coef*z3);V33_z4  <-  -1*sum(V33_comp.coef*z4)
  V33_z5  <-  -1*sum(V33_comp.coef*z5);V33_z6  <-  -1*sum(V33_comp.coef*z6);V33_z7  <-  -1*sum(V33_comp.coef*z7)
  
  V33_z11  <-  -1*sum(V33_comp.coef*z1*z1) ##V_z11 = V_z1z1
  V33_z12  <-  -1*sum(V33_comp.coef*z1*z2);V33_z13  <-  -1*sum(V33_comp.coef*z1*z3);V33_z14  <-  -1*sum(V33_comp.coef*z1*z4)
  V33_z15  <-  -1*sum(V33_comp.coef*z1*z5);V33_z16  <-  -1*sum(V33_comp.coef*z1*z6);V33_z17  <-  -1*sum(V33_comp.coef*z1*z7)
  
  V33_z22  <-  -1*sum(V33_comp.coef*z2*z2) ##V_z22 = V_z2z2
  V33_z23  <-  -1*sum(V33_comp.coef*z2*z3);V33_z24  <-  -1*sum(V33_comp.coef*z2*z4)
  V33_z25  <-  -1*sum(V33_comp.coef*z2*z5);V33_z26  <-  -1*sum(V33_comp.coef*z2*z6);V33_z27  <-  -1*sum(V33_comp.coef*z2*z7)
  
  V33_z33  <-  -1*sum(V33_comp.coef*z3*z3);V33_z34  <-  -1*sum(V33_comp.coef*z3*z4)
  V33_z35  <-  -1*sum(V33_comp.coef*z3*z5);V33_z36  <-  -1*sum(V33_comp.coef*z3*z6);V33_z37  <-  -1*sum(V33_comp.coef*z3*z7)
  
  
  V33_z44  <-  -1*sum(V33_comp.coef*z4*z4)
  V33_z45  <-  -1*sum(V33_comp.coef*z4*z5);V33_z46  <-  -1*sum(V33_comp.coef*z4*z6);V33_z47  <-  -1*sum(V33_comp.coef*z4*z7)
  
  V33_z55  <-  -1*sum(V33_comp.coef*z5*z5);V33_z56  <-  -1*sum(V33_comp.coef*z5*z6);V33_z57  <-  -1*sum(V33_comp.coef*z5*z7)
  
  V33_z66  <-  -1*sum(V33_comp.coef*z6*z6);V33_z67  <-  -1*sum(V33_comp.coef*z6*z7)
  
  V33_z77  <-  -1*sum(V33_comp.coef*z7*z7)
  
  til_V33 <- matrix(NA,8,8)
  top_tri <- matrix(c(V33_11,V33_z1,V33_z2,V33_z3,V33_z4,V33_z5,V33_z6,V33_z7,
                      NA,V33_z11,V33_z12,V33_z13,V33_z14,V33_z15,V33_z16,V33_z17,
                      NA,NA,V33_z22,V33_z23,V33_z24,V33_z25,V33_z26,V33_z27,
                      NA,NA,NA,V33_z33,V33_z34,V33_z35,V33_z36,V33_z37,
                      NA,NA,NA,NA,V33_z44, V33_z45, V33_z46, V33_z47,
                      NA,NA,NA,NA,NA,V33_z55,V33_z56,V33_z57,
                      NA,NA,NA,NA,NA,NA,V33_z66,V33_z67,
                      NA,NA,NA,NA,NA,NA,NA, V33_z77),8,8,byrow = TRUE)
  for(i in 1:8){
    for(j in 1:8){
      if(i>=j){
        til_V33[i,j]<- top_tri[j,i]
      }else{
        til_V33[i,j]<- top_tri[i,j]
      }
    }
  }
  
  til_V33.invers <- solve(til_V33)
  til_phi <- sum(til_pi/til_fbig)
  
  N_asyVar <- c(til_phi - 1 - t(as.matrix(til_V23)) %*%  til_V33.invers %*% as.matrix(til_V23))
  alpha_asyVar <- til_alpha*(1-til_alpha)
  beta_asyVar <-  -1* til_V33.invers
  
  
  ## -----------------------------   overall&grouped pointwise confidence band ------------------------------------#
  Z_Sig <- qnorm(1-Siglev/2)  
  Delta2 <- function(t){
   
    fbig_t <- NULL
    fbig1_t <- NULL
    ter_coef <- NULL
    for (i in 1:n) {
      qz_i <- c(1,z1[i],z2[i],z3[i],z4[i],z5[i],z6[i],z7[i])
      
      fbig_t[i] <- f_b(t,til_theta[i])
      fbig1_t[i] <- f1_b(t,til_theta[i])
      
      ter_coef <- rbind(ter_coef,
                        (fbig1_t[i] - fbig_t[i]* til_f1.big[i]/til_fbig[i]  )*til_theta[i]*til_pi[i]*qz_i )
      
    }
   apply(ter_coef, 2, sum)
  }

  
  #------------------------- overall&grouped distribution  ----------------------------#
  idx_size1 <- (z3 ==1 &  z4 ==0)
  idx_size2 <- (z3 ==0 &  z4 ==1)
  idx_size3 <- (z3 ==0 &  z4 ==0)
  
  fz1_f.til <- til_pi*idx_size1/sum( til_pi*idx_size1 ) ## size = "<20"density
  fz2_f.til <- til_pi*idx_size2/sum( til_pi*idx_size2 )
  fz3_f.til <- til_pi*idx_size3/sum( til_pi*idx_size3 )
  
  
  fb.band <- function(t){
    fb_val <- NULL
    for (k in 1:n) {
      fb_val[k] <- f_b(t,til_theta[k])
    }
    f_big        <- sum(til_pi*fb_val )
    fb_val.size1 <- sum(fz1_f.til   * fb_val )
    fb_val.size2 <- sum(fz2_f.til   * fb_val )
    fb_val.size3 <- sum(fz3_f.til   * fb_val )
    
    critcal <- sqrt(  t(as.matrix(Delta2(t)))%*% beta_asyVar %*% as.matrix(Delta2(t)) ) /sqrt(til_N)
    
    c(f_big, fb_val.size1, fb_val.size2,fb_val.size3, critcal)
    
  }
  
  rfst <- seq(0,20,by = 0.1)
  ter_val    <- NULL
  for (i in 1: length(rfst)) {
    ter_val <- rbind(ter_val,  fb.band(rfst[i]) )
  }
  
  fbig_int.or <- data.frame(fbig = ter_val[,1], fbig.u = ter_val[,1] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,1] -  Z_Sig* ter_val[,5])  ## overall pointwise confidence band 
  
  fbig_int.s1 <- data.frame(fbig = ter_val[,2], fbig.u = ter_val[,2] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,2] -  Z_Sig* ter_val[,5])  ## Size1 
  
  fbig_int.s2 <- data.frame(fbig = ter_val[,3], fbig.u = ter_val[,3] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,3] -  Z_Sig* ter_val[,5])  ## Size2
  fbig_int.s3 <- data.frame(fbig = ter_val[,4], fbig.u = ter_val[,4] +  Z_Sig* ter_val[,5] ,
                            fbig.l= ter_val[,4] -  Z_Sig* ter_val[,5])  ## Size3
  
  
  ##-------------------output -------------------------#
  list("N_asyVar"= N_asyVar, "alpha_asyVar" = alpha_asyVar, "beta_asyVar"=beta_asyVar,"fbig_intval" =  fbig_int.or,
           "fbig_intval.size1" = fbig_int.s1,"fbig_intval.size2" = fbig_int.s2,"fbig_intval.size3" = fbig_int.s3  )
}



#----------------------------------------------------------------------------------------#

#       find the margin cumulative distribution function of rfstime                       #                      

#----------------------------------------------------------------------------------------#

fb_f <- function(t,OBS=obs,intib){
  n <- nrow(OBS)
  x <- OBS[,1]
  y <- OBS[,2]
  z1 <- OBS[,3]
  z2 <- OBS[,4]
  z3 <- OBS[,5]
  z4 <- OBS[,6]
  z5 <- OBS[,7]
  z6 <- OBS[,8]
  z7 <- obs[,9]
  inti0 = intib
  
  
  hat_out <- maxlike.par(OBS,inti = inti0) # optimal intinal value : rep(0.01,5)
  hat_N <- hat_out[[1]][1]
  hat_alpha <- hat_out[[1]][2]
  hat_beta <- hat_out[[1]][3:10]
  hat_lab <- hat_out[[1]][11]
  
  com <- NULL
  hat_theta <- NULL
  hat_pi <- NULL
  for (i in 1:n) {
    hat_theta[i] <- exp(hat_beta[1] + hat_beta[2]*z1[i]+ hat_beta[3]*z2[i]+ hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                        + hat_beta[6]*z5[i]+ hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    hat_pi[i] <- 1/(1  + hat_lab*(f_b(y[i],hat_theta[i]) - hat_alpha))/n
    com <- rbind(com,  hat_pi[i]*f_b(t,hat_theta[i]) )
  }
  sum(com)
}


fb_c <- function(t,OBS=obs,intib){
  n <- nrow(OBS)
  x <- OBS[,1]
  y <- OBS[,2]
  z1 <- OBS[,3]
  z2 <- OBS[,4]
  z3 <- OBS[,5]
  z4 <- OBS[,6]
  z5 <- OBS[,7]
  z6 <- OBS[,8]
  z7 <- obs[,9]
  
  inti0 = intib
  
  til_out <- max.clike(OBS,inti = inti0) 
  til_N <- til_out[[1]][1]
  til_alpha <- til_out[[1]][2]
  til_beta  <- til_out[[1]][3:10]
  
  com <- NULL
  til_theta <- NULL
  til_pi <- NULL
  for (i in 1:n) {
    til_theta[i] <- exp(til_beta[1] + til_beta[2]*z1[i]+ til_beta[3]*z2[i]+ til_beta[4]*z3[i]+ til_beta[5]*z4[i]
                        + til_beta[6]*z5[i]+ til_beta[7]*z6[i]+ til_beta[8]*z7[i])
    til_pi[i] <- 1/f_b(y[i],til_theta[i])/til_N
    com <- rbind(com, til_pi[i]*f_b(t,til_theta[i]) )
  }
  sum(com)
}




##------------------------------ group by size  ----------------------------------------------##
#
#                           size = "<20":   size1=1,size2=0
#                           size = "20-50": size1=0,size2=1
#                           size = ">50" :  size1=0,size2=0
#
#------------------------------------------------------------------------------------------------------##

fb_f.size123 <- function(t,OBS = obs,intib){
  n <- nrow(OBS)
  x <- OBS[,1]
  y <- OBS[,2]
  z1 <- OBS[,3]
  z2 <- OBS[,4]
  z3 <- OBS[,5]##size1 = "<20"
  z4 <- OBS[,6]##size2 = "20-50"
  z5 <- OBS[,7]
  z6 <- OBS[,8]
  z7 <- OBS[,9]
  
  
  inti0 = intib
  
  idx_size1 <- (z3 ==1 &  z4 ==0)
  idx_size2 <- (z3 ==0 &  z4 ==1)
  idx_size3 <- (z3 ==0 &  z4 ==0)
  
  hat_out <- maxlike.par(obs,inti = inti0) # optimal intinal value : rep(0.01,5)
  hat_N <- hat_out[[1]][1]
  hat_alpha <- hat_out[[1]][2]
  hat_beta <- hat_out[[1]][3:10]
  hat_lab <- hat_out[[1]][11]
  
  fb_f.hat <- NULL
  hat_theta <- NULL
  hat_pi <- NULL
  for (i in 1:n) {
    hat_theta[i] <- exp(hat_beta[1] + hat_beta[2]*z1[i]+ hat_beta[3]*z2[i]+ hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                        + hat_beta[6]*z5[i]+ hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    hat_pi[i] <- 1/(1  + hat_lab*(f_b(y[i],hat_theta[i]) - hat_alpha))/n
    fb_f.hat[i] <- f_b(t,hat_theta[i])
  }
  
  fz1_f.hat <- hat_pi*idx_size1/sum( hat_pi*idx_size1 )
  fz2_f.hat <- hat_pi*idx_size2/sum( hat_pi*idx_size2 )
  fz3_f.hat <- hat_pi*idx_size3/sum( hat_pi*idx_size3 )
  
  
  
  fb_f.size1 <- sum(fz1_f.hat* fb_f.hat )
  fb_f.size2 <- sum(fz2_f.hat* fb_f.hat )
  fb_f.size3 <- sum(fz3_f.hat* fb_f.hat )
  
  c(fb_f.size1,fb_f.size2,fb_f.size3)
  
}


fb_c.size123 <- function(t,OBS = obs,intib){
  n <- nrow(OBS)
  x <- OBS[,1]
  y <- OBS[,2]
  z1 <- OBS[,3]
  z2 <- OBS[,4]
  z3 <- OBS[,5]##size1 = "<20"
  z4 <- OBS[,6]##size2 = "20-50"
  z5 <- OBS[,7]
  z6 <- OBS[,8]
  z7 <- OBS[,9]
  
  inti0 = intib
  
  idx_size1 <- (z3 ==1 &  z4 ==0)
  idx_size2 <- (z3 ==0 &  z4 ==1)
  idx_size3 <- (z3 ==0 &  z4 ==0)
  
  til_out <- max.clike(obs,inti = inti0) 
  til_N <- til_out[[1]][1]
  til_alpha <- til_out[[1]][2]
  til_beta  <- til_out[[1]][3:10]
  
  
  fb_c.til <- NULL
  til_theta <- NULL
  til_pi <- NULL
  for (i in 1:n) {
    til_theta[i] <- exp(til_beta[1] + til_beta[2]*z1[i]+ til_beta[3]*z2[i]+ til_beta[4]*z3[i]+ til_beta[5]*z4[i]
                        + til_beta[6]*z5[i]+ til_beta[7]*z6[i]+ til_beta[8]*z7[i])
    til_pi[i] <- 1/f_b(y[i],til_theta[i])/til_N
    fb_c.til[i] <- f_b(t, til_theta[i])
  }
  
  fz1_f.til <- til_pi*idx_size1/sum( til_pi*idx_size1 )
  fz2_f.til <- til_pi*idx_size2/sum( til_pi*idx_size2 )
  fz3_f.til <- til_pi*idx_size3/sum( til_pi*idx_size3 )
  
  
  fb_c.size1 <- sum(fz1_f.til* fb_c.til )
  fb_c.size2 <- sum(fz2_f.til* fb_c.til )
  fb_c.size3 <- sum(fz3_f.til* fb_c.til )
  
  c(fb_c.size1,fb_c.size2,fb_c.size3)
}