
#---------------------   call for real_Main.R   ---------------------------------------------#
#
#     data: rotterdam breast cancer data
#           Only patients with positive lymph nodes were considered for rfstime distribution
#          rfstime refers to the time from the surgical procedure to the first recurrence or immediate death without recurrence        
#  
#------------------------------------------------------------------------------------------#

library(survival)

brecan0  <- rotterdam
brecan0$rfstime <- pmin(rotterdam$rtime, rotterdam$dtime)/365.25 #years
brecan0$status  <- pmax(rotterdam$recur, rotterdam$death)        #1: failure, failure means dead or recurrent
brecan0$mxfyea  <- pmax(rotterdam$rtime, rotterdam$dtime)/365.25 #Maximum follow-up year
brecan0$maxftoyear <- brecan0$year + brecan0$mxfyea; range(brecan0$maxftoyear) ##Year of maximum follow-up
brecan0$size1 <-   (rotterdam$size == "<=20")+0 ## size > 50 is baseline
brecan0$size2 <-   (rotterdam$size == "20-50")+0


indicator <- (brecan0$nodes > 0)
brecan <- brecan0[indicator,]
N <- nrow(brecan);











##- ------------------------------------------------------------------------------------------#
#                           right truncation                                                                  
#                                                             
##------------------------------------------------------------------------------------------- #
start <- 1978
end   <- 1993.5
trunct <- end - brecan$year
idx  <- ( brecan$year >= 1978 &  brecan$rfstime <= trunct & brecan$status ==1)

brecan_obs <- brecan[idx,]
alpha_t    <- sum(idx)/N; alpha_t 

x <- brecan_obs$rfstime
y <- end - brecan_obs$year
z1 <- brecan_obs$chemo 
z2 <- brecan_obs$hormon
z3 <- brecan_obs$size1
z4 <- brecan_obs$size2
z5 <- (brecan_obs$age > 50)+0  
z6 <- log(1+brecan_obs$nodes)
z7 <- log(1+sqrt(brecan_obs$er*brecan_obs$age))

obs0 <- cbind(x,y,z1,z2,z3,z4,z5,z6,z7)
obs <- obs0[order(obs0[,1]),] #Sort the dataset in ascending order by rfstime
n <- nrow(obs)

##-------------------------  Find the best initial value   ------------------------------#
maxlike_f <- NULL
maxlike_c <- NULL
inti_f <- NULL
inti_c <- NULL
par_fest <- NULL
par_cest <- NULL
for (i in 1:100) {

  intif <- runif(8,-1,1)
  intic <- runif(8,-1,1)
  
  
  resf <- maxlike.par(obs,inti = intif)
  resc <- max.clike(obs, inti = intic)
  
  par_fest <- rbind(par_fest,resf[[1]])
  par_cest <- rbind(par_cest,resc[[1]])
  
  maxlike_f <- rbind(maxlike_f,resf[[2]])
  maxlike_c <- rbind(maxlike_c,resc[[2]])
  
  inti_f <- rbind(inti_f,intif)
  inti_c <- rbind(inti_c,intic)
  print(i)
}


idxf <- which(maxlike_f == max(maxlike_f) -1 ) 
idxc <- which(maxlike_c == max(maxlike_c) -1 )
options(digits = 10)
intif <- inti_f[idxf,];intif
intic <- inti_c[idxc,];intic

par_fest[idxf,] 
par_cest[idxc,] 
maxlike_f[idxf]
maxlike_c[idxc]





#---------------------------------------------------------------------------------#
#  Model test :   call for Final_Main.R  &  Final_Test.R
#            hat_RS = 1.58219
#
#           qchibarsq(1-0.05, df=1) = 2.705543
#----------------------------------------------------------------------------------#
intif_best <- c(-0.58622332476, -0.92664983822, -0.49120514188, -0.30701229675,  0.05260184128, -0.19268889073,0.54528445425, -0.44097818621)#1993.5
intic_best <- c(-0.6681207479,  -0.7034504893,  -0.5959785758 , -0.9838226857 ,  0.6489063473 , -0.3959745718 , 0.1593918940, -0.2778926282)#1993.5

hat_par <- maxlike.par(obs = obs,inti = intif_best)[[1]];hat_par
til_par <- max.clike(obs = obs,inti = intic_best)[[1]];til_par

hat_cb     <- cb_est(obs,hat_par )
hat_inV33  <- inV33_est(obs,hat_par )
hat_uts    <- U_ts(obs,hat_par)
hat_RS  <- max(hat_uts,0)^2;hat_RS

c( (hat_RS > qchibarsq(1-0.01,df=1)), (hat_RS > qchibarsq(1-0.05, df = 1))  ) 


#-------------------------------------------------------------------------------------#
#  Model :  call for  Final_Test.R
#                
# --------------------------------------------------------------------------------------#


asy_res.f <- Asyvar.full(obs = obs,hat_par = hat_par ,0.05)
asy_res.c <- Asyvar.con(obs = obs,til_par =til_par,0.05)

sqrt(asy_res.f$N_asyVar);  sqrt(asy_res.c$N_asyVar)
sqrt(asy_res.f$alpha_asyVar);  sqrt(asy_res.c$alpha_asyVar)
sqrt( diag( asy_res.f$beta_asyVar )) ; sqrt( diag(asy_res.c$beta_asyVar))



##--------------------  95% confidence band    ----------------------------------¡ª¡ª# 

hat_N <- hat_par[1]
hat_Nsd <- sqrt(asy_res.f$N_asyVar)

hat_N - qnorm(1-0.05/2)*sqrt(hat_N )* hat_Nsd 
hat_N + qnorm(1-0.05/2)*sqrt(hat_N )* hat_Nsd 


til_N <- til_par[1]
til_Nsd <- sqrt(asy_res.c$N_asyVar)

til_N - qnorm(1-0.05/2)*sqrt(til_N )* til_Nsd 
til_N + qnorm(1-0.05/2)*sqrt(til_N )* til_Nsd 


hat_alpha <- hat_par[2]
hat_Asd <- sqrt(asy_res.f$alpha_asyVar)

hat_alpha - qnorm(1-0.05/2)*hat_Asd/sqrt(hat_N)
hat_alpha + qnorm(1-0.05/2)*hat_Asd/sqrt(hat_N)


til_alpha <- til_par[2]
til_Asd <- sqrt(asy_res.c$alpha_asyVar)

til_alpha - qnorm(1-0.05/2)*til_Asd/sqrt(til_N)
til_alpha + qnorm(1-0.05/2)*til_Asd/sqrt(til_N)

