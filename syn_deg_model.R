# Model of synthesis and degredation in the model ------------------------------

time_final <- 120
initial <- 1000

deg_k <- 0.05
syn <- 50

time <- c(0:time_final)

# Old model style

dat <- NULL

for (i in 1:length(time)){
  
  temp <- initial / deg_k
  dat$temp[i] <- temp / deg_k
  
}

dat <- as.data.frame(dat)



# One phase decay Y = (Y0 - Plateau) * e(-K * X) + Plateau



  



