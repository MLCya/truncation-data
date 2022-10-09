# ----------------------        call for  Final_Main.R ----------------------------------------#


#-------------------------------------------------------------------------#

#    generate the quantiles of  "mixed" chi-squared distributions         #
#        (mixtures of chi-square(n) and chi-square(n-1)                   #

#-------------------------------------------------------------------------#
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





#----------------------------------------------------------------------------#

#  theta = exp(beta %*% q(z))    q(z) = (1,z1,z2,z3,z4,z5,z6,z7) 
#  target variable x:rfstime  follows exponential distribution
#----------------------------------------------------------------------------#

cb_est <- function(obs,hat_par){
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
  hat_com <- NULL
  for (i in 1:n) {
    qz_i <- c(1, z1[i], z2[i],z3[i],z4[i],z5[i],z6[i],z7[i])
    hat_theta[i] <-  exp(hat_beta[1]+hat_beta[2]*z1[i] + hat_beta[3]*z2[i] + hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                         + hat_beta[6]*z5[i] + hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    hat_pi <- 1/n * 1/(1+hat_lab*(f_b(y[i], hat_theta[i]) - hat_alpha )   )
    fn <- function(x){
      (f2_s(x,hat_theta[i])*hat_theta[i]+f1_s(x,hat_theta[i]))*(1/hat_theta[i] - x)
    }
    
    int_val <- integrate(fn,lower = 0,upper = y[i])[[1]]
    sec_val <- (f2_b(y[i],hat_theta[i])*hat_theta[i]+f1_b(y[i],hat_theta[i]))*f1_b(y[i],hat_theta[i])/f_b(y[i],hat_theta[i])
    
    hat_com <- rbind(hat_com,  (int_val - sec_val)*hat_theta[i]^2*hat_pi*qz_i)
  }
  -1*apply(hat_com, 2, sum)
}

#--------------------------------------------------------------------------------------#

##      estimate V_{33}^{-1} when q(z) = (1,z1,z2,z3,z4,z5,z6,z7)                         #
##                                               #

#-------------------------------------------------------------------------------------#
inV33_est <- function(obs,hat_par){
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
  
  hat_theta <- NULL  ## h(z,til_beta)
  hat_coef <- NULL
  for (i in 1:n) {
    hat_theta[i] <-  exp(hat_beta[1]+hat_beta[2]*z1[i] + hat_beta[3]*z2[i] + hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                         + hat_beta[6]*z5[i] + hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    hat_pi <- 1/n * 1/(1+hat_lab*(f_b(y[i], hat_theta[i]) - hat_alpha )   )
    
    fn <- function(x){
      (1/hat_theta[i] - 2*x +hat_theta[i]*x^2)/exp(hat_theta[i]*x)
    }
    int_val <- integrate(fn,lower = 0,upper = y[i])[[1]]
    sec_val <- f1_b(y[i],hat_theta[i])^2/f_b(y[i], hat_theta[i])
    hat_coef[i] <- (int_val - sec_val )*hat_theta[i]^2*hat_pi
  }  
  
  
  ###       The  components of V_{33}        ##########
  V33_11  <-  -1*sum(hat_coef*1*1)
  
  V33_z1  <-  -1*sum(hat_coef*z1) ##V_z1
  V33_z2  <-  -1*sum(hat_coef*z2);V33_z3  <-  -1*sum(hat_coef*z3);V33_z4  <-  -1*sum(hat_coef*z4)
  V33_z5  <-  -1*sum(hat_coef*z5);V33_z6  <-  -1*sum(hat_coef*z6);V33_z7  <-  -1*sum(hat_coef*z7)
  
  V33_z11  <-  -1*sum(hat_coef*z1*z1) ##V_z11 = V_z1z1
  V33_z12  <-  -1*sum(hat_coef*z1*z2);V33_z13  <-  -1*sum(hat_coef*z1*z3);V33_z14  <-  -1*sum(hat_coef*z1*z4)
  V33_z15  <-  -1*sum(hat_coef*z1*z5);V33_z16  <-  -1*sum(hat_coef*z1*z6);V33_z17  <-  -1*sum(hat_coef*z1*z7)
  
  V33_z22  <-  -1*sum(hat_coef*z2*z2) ##V_z22 = V_z2z2
  V33_z23  <-  -1*sum(hat_coef*z2*z3);V33_z24  <-  -1*sum(hat_coef*z2*z4)
  V33_z25  <-  -1*sum(hat_coef*z2*z5);V33_z26  <-  -1*sum(hat_coef*z2*z6);V33_z27  <-  -1*sum(hat_coef*z2*z7)
  
  V33_z33  <-  -1*sum(hat_coef*z3*z3);V33_z34  <-  -1*sum(hat_coef*z3*z4)
  V33_z35  <-  -1*sum(hat_coef*z3*z5);V33_z36  <-  -1*sum(hat_coef*z3*z6);V33_z37  <-  -1*sum(hat_coef*z3*z7)
  
  
  V33_z44  <-  -1*sum(hat_coef*z4*z4)
  V33_z45  <-  -1*sum(hat_coef*z4*z5);V33_z46  <-  -1*sum(hat_coef*z4*z6);V33_z47  <-  -1*sum(hat_coef*z4*z7)
  
  V33_z55  <-  -1*sum(hat_coef*z5*z5);V33_z56  <-  -1*sum(hat_coef*z5*z6);V33_z57  <-  -1*sum(hat_coef*z5*z7)
  
  V33_z66  <-  -1*sum(hat_coef*z6*z6);V33_z67  <-  -1*sum(hat_coef*z6*z7)
  
  V33_z77  <-  -1*sum(hat_coef*z7*z7)
  
  V33_est <- matrix(NA,8,8)
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
        V33_est[i,j]<- top_tri[j,i]
      }else{
        V33_est[i,j]<- top_tri[i,j]
      }
    }
  }
  solve(V33_est)
}


#--------------------------------------------------------------------------------------#

##      estimate Sigam^2, when q(z) = (1,z1,z2,z3,z4,z5,z6,z7)      #

#-------------------------------------------------------------------------------------#

############################ to estimate sigam_s^2 #############################
Sigmas2_est <- function(obs,hat_par){
  n = nrow(obs)
  x  <- obs[,1]
  y  <- obs[,2]
  z1 <- obs[,3]
  z2 <- obs[,4]
  z3 <- obs[,5]##size1 = "<20"
  z4 <- obs[,6]##size2 = "20-50"
  z5 <- obs[,7]
  z6 <- obs[,8]
  z7 <- obs[,9]
  
  
  hat_cb <- cb_est(obs, hat_par)
  hat_inV33 <- inV33_est(obs, hat_par )
  cv   <- t(as.matrix(hat_cb)) %*% hat_inV33
  
  
  hat_N   <- hat_par[1]
  hat_alpha <- hat_par[2]
  hat_beta <- matrix(hat_par[3:10],ncol = 8,dimnames = NULL)
  hat_lab  <-  hat_par[11]
  
  
  hat_theta <- NULL  ## h(z,til_beta)
  hat_coef <- NULL
  for (i in 1:n) {
    qz_i <- c(1,z1[i],z2[i],z3[i],z4[i],z5[i],z6[i],z7[i])
    hat_theta[i] <-  exp(hat_beta[1]+hat_beta[2]*z1[i] + hat_beta[3]*z2[i] + hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                         + hat_beta[6]*z5[i] + hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    hat_pi <- 1/n * 1/(1+hat_lab*(f_b(y[i], hat_theta[i]) - hat_alpha )   )
    
    cvq <- sum( cv*qz_i)
    fn <- function(x){
      ((hat_theta[i]*x-2)*hat_theta[i]*x+(1-hat_theta[i]*x)*(1-cvq))^2/(hat_theta[i]*exp(hat_theta[i]*x))
    }
    
    int_val <- integrate(fn,lower = 0,upper = y[i])[[1]]
    sec_val <- (f2_b(y[i],hat_theta[i])*hat_theta[i]+f1_b(y[i],hat_theta[i])*(1-cvq))^2/f_b(y[i],hat_theta[i]) 
    hat_coef[i] <- (int_val - sec_val )*hat_theta[i]^2 * hat_pi
  }
  sum(hat_coef)
}

######################     the standard score statistics     ###################################
U_ts <- function(obs,hat_par){
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
  
  hat_sigms <- sqrt(Sigmas2_est(obs,hat_par )) ##  standard deviation estimator
  
  hat_theta <- NULL  ## h(z,til_beta)
  terp <- NULL
  for (i in 1:n) {
    hat_theta[i] <-  exp(hat_beta[1]+hat_beta[2]*z1[i] + hat_beta[3]*z2[i] + hat_beta[4]*z3[i]+ hat_beta[5]*z4[i]
                         + hat_beta[6]*z5[i] + hat_beta[7]*z6[i]+ hat_beta[8]*z7[i])
    #fir_val <- hat_theta[i]*x[i]*(hat_theta[i]*x[i] -2) + (1- hat_theta[i]*x[i])
    #sec_val <- hat_theta[i]*y[i]*(1 - hat_theta[i]*y[i])/(exp(hat_theta[i]*y[i])-1)
    #terp[i] <- fir_val - sec_val
    fir_val <- (f2_s(x[i],hat_theta[i])*hat_theta[i] + f1_s(x[i],hat_theta[i]))/f_s(x[i], hat_theta[i])
    sec_val <- (f2_b(y[i],hat_theta[i])*hat_theta[i] + f1_b(y[i],hat_theta[i]))/f_b(y[i],hat_theta[i])
    terp[i] <- (fir_val - sec_val)*hat_theta[i]
  }
  1/sqrt(hat_N)*sum(terp)/hat_sigms
}