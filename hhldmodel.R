####################################################################################################
# Build a simple model to describe yield responses to variations in pipe irrigation 
# Author: Drew Gower
####################################################################################################
## Clear Console and load required packages
rm(list=ls(all=TRUE))
library(abind)


## Start time counter
ptm <- proc.time()


## Define Parameter values 
# Climate Parameters
alpha_mn <- 0.015 # m
alpha_cv <- 0.25 # dim
lambda_mn <- 0.3 # day^-1
lambda_cv <- 0.25 # dim
Emax <- 0.005 # m/day
Ew <- 0.0001 # m/day
Tseas <- 110 # days

# Household Parameters
frac_sub.v <- seq(0, 1, 0.2) # dim
A_tot <- 5000 # m^2
Ymin <- 1000 # kg

# Crop Parameters
Zr <- 0.5 # m
Ymax_sub <- 0.3 # kg/m^2
Ymax_mar <- 0.3 # kg/m^2
q_sub <- 2 # dim
q_mar <- q_sub / 2 # dim
k_sub <- 1 # dim 
k_mar <- k_sub / 2 # m 
r_sub <- 0.5 # dim
r_mar <- 0.5 # dim

# Soil Prameters (Sandy Loam)
Ks <- 0.8 # m/day
n <- 0.43 # dim
beta <- 13.8 # dim
sh <- 0.14 # dim
sw <- 0.18 # dim
ss <- 0.46 # dim
sfc <- 0.56 # dim

# Economic Parameters
cy_sub <- 0.2 # $/kg
cy_mar <- cy_sub * 1.4 # $/kg
cw <- 0.02 # $/day

# Irrigation Parameters
smin <- ss # dim
smax <- sfc # dim
Qrte <- 10 # m^3/day
var.v <- c(1/1000, 1/3, 1)

## Begin Simulation
trials <- 10000
buffer <- 10 # days
interval <- 24 # intervals/day

# Define Superduperglobal Data Structures
P_tot.m <- matrix(NA,interval*(Tseas+buffer),trials)
S_tot.a1 <- array(NA,c(interval*Tseas,trials,length(frac_sub.v),length(var.v)))
S_tot_dens.l1 <- vector("list",length(var.v))
S_tot_dens.a <- array(NA,c(3,length(frac_sub.v),length(var.v)))
I_tot.a1 <- array(NA,c(Tseas,trials,length(frac_sub.v),length(var.v)))
I_tot_dens.l1 <- vector("list",length(var.v))
I_tot_dens.a <- array(NA,c(3,length(frac_sub.v),length(var.v)))
Y_sub.a <- array(NA,c(trials,length(frac_sub.v),length(var.v)))
Y_mar.a <- array(NA,c(trials,length(frac_sub.v),length(var.v)))
Y_sub_dens.l1 <- vector("list",length(var.v))
Y_mar_dens.l1 <- vector("list",length(var.v))
Y_com_dens.a <- array(NA,c(3,length(frac_sub.v),length(var.v)))
Y_tot.a <- array(NA,c(trials,length(frac_sub.v),length(var.v)))
Y_tot_dens.l1 <- vector("list",length(var.v))
Y_tot_dens.a <- array(NA,c(3,length(frac_sub.v),length(var.v)))
G_tot.a <- array(NA,c(trials,length(frac_sub.v),length(var.v)))
G_tot_dens.l1 <- vector("list",length(var.v))
G_tot_dens.a <- array(NA,c(3,length(frac_sub.v),length(var.v)))
stat_tot.a <- stat_10.a <- stat_50.a <- stat_90.a <- 
  array(NA,c(2,length(frac_sub.v),length(var.v)))

# Define Rainfall Pattern
for(i in 1:trials){
  
  # Pull Seasonal Alpha and Lambda from Gamma Distribution 
  alpha <- rgamma(1,shape=1/alpha_cv^2,scale=alpha_mn*alpha_cv^2)
  lambda <- rgamma(1,shape=1/lambda_cv^2,scale=lambda_mn*lambda_cv^2)
  
  # Define Rainfall for each Interval
  P_days.v <- rbinom(interval*(Tseas+buffer),1,lambda/interval)
  P_amt.v <- rexp(length(P_days.v[P_days.v>0]),1/alpha)
  P_tot.v <- replace(P_days.v,P_days.v==1,P_amt.v)
  P_tot.v[P_tot.v<0] <- 0
  P_tot.m[,i] <- P_tot.v
}

# Begin First iteration
for(g in 1:length(var.v)) {
  
  Qfre <- var.v[g]
  
  # Define Superglobal Data Structures
  S_tot.a2 <- array(NA,c(interval*Tseas,trials,length(frac_sub.v)))
  S_tot_dens.l2 <- vector("list",length(frac_sub.v))
  S_tot_dens.m <- matrix(NA,3,length(frac_sub.v))
  I_tot.a2 <- array(NA,c(Tseas,trials,length(frac_sub.v)))
  I_tot_dens.l2 <- vector("list",length(frac_sub.v))
  I_tot_dens.m <- matrix(NA,3,length(frac_sub.v))
  Y_sub.m <- matrix(NA,trials,length(frac_sub.v))
  Y_mar.m <- matrix(NA,trials,length(frac_sub.v))
  Y_sub_dens.l2 <- vector("list",length(frac_sub.v))
  Y_mar_dens.l2 <- vector("list",length(frac_sub.v))
  Y_com_dens.m <- matrix(NA,3,length(frac_sub.v))
  Y_tot.m <- matrix(NA,trials,length(frac_sub.v))
  Y_tot_dens.l2 <- vector("list",length(frac_sub.v))
  Y_tot_dens.m <- matrix(NA,3,length(frac_sub.v))
  G_tot.m <- matrix(NA,trials,length(frac_sub.v))
  G_tot_dens.l2 <- vector("list",length(frac_sub.v))
  G_tot_dens.m <- matrix(NA,3,length(frac_sub.v))
  stat_tot.m <- matrix(NA,2,length(frac_sub.v))
  
  for(h in 1:length(frac_sub.v)) {
    
    # Complete Area Calculations
    A_sub <- frac_sub.v[h]*A_tot # m^2
    A_mar <- A_tot-A_sub # m^2
    
    # Define Global Data Structures
    S_tot.m <- matrix(NA,interval*Tseas,trials)
    I_tot.m <- matrix(NA,Tseas,trials)
    ET_tot.m <- matrix(NA,Tseas,trials)
    R_tot.m <- matrix(NA,Tseas,trials)
    L_tot.m <- matrix(NA,Tseas,trials)
    theta_sub.v <- vector("numeric",trials)
    theta_mar.v <- vector("numeric",trials)
    Y_sub.v <- vector("numeric",trials)
    Y_mar.v <- vector("numeric",trials)
    Y_tot.v <- vector("numeric",trials)
    G_tot.v <- vector("numeric",trials)
    
    # Loop through trials
    for(i in 1:trials) {
      
      # Create Daily and Subdaily Data Storage
      P_tot.v <- P_tot.m[,i]
      S_tot.v <- vector("numeric",interval*(Tseas+buffer))
      I_tot.v <- vector("numeric",Tseas+buffer)
      ET_tot.v <- vector("numeric",Tseas+buffer)
      R_tot.v <- vector("numeric",Tseas+buffer)
      L_tot.v <- vector("numeric",Tseas+buffer)
      c_day.v <- vector("numeric",Tseas+buffer)
      
      # Define Initial Values
      S_tot.old <- sw # dim
      
      # Loop through days
      cnt <- 0
      for(j in 1:(buffer+Tseas)) {
        S_tot <- S_tot.old
        
        # Create Temporary Data Storage
        S_tot.temp <- vector("numeric",interval)
        I_tot.temp <- vector("numeric",interval)
        ET_tot.temp <- vector("numeric",interval)
        R_tot.temp <- vector("numeric",interval)
        L_tot.temp <- vector("numeric",interval)
        
        # Determine if irrigation is available
        if(j%%(1/Qfre)==0) {
          c_day <- cw
          Q_day <- Qrte
        } else {
          c_day <- 0
          Q_day <- 0  
        }
        
        # Begin Subdaily Iteration
        for(k in 1:interval) {
          
          # Add Rainfall
          S_tot <- S_tot+P_tot.v[(j-1)*interval+k]/(n*Zr)
          
          # Remove Surface Runoff
          if(S_tot>1) {
            R_tot <- S_tot-1
            S_tot <- 1
          } else {
            R_tot <- 0
          }
          
          # Remove Leakage
          if(S_tot>sfc) {
            L_tot <- (Ks/interval)*(exp(beta*(S_tot-sfc))-1)/(exp(beta*(1-sfc))-1)
            S_tot <- S_tot-L_tot/(n*Zr)
          } else {
            L_tot <- 0
          }
          
          # Remove Evapotranspiration Losses
          if(S_tot>ss) {
            ET_tot <- Emax/interval
          } else if(S_tot>sw && S_tot<=ss) {
            ET_tot <- (Ew/interval)+((Emax-Ew)/interval)*(S_tot-sw)/(ss-sw)
          } else if(S_tot>sh && S_tot<=sw) {
            ET_tot <- (Ew/interval)*(S_tot-sh)/(sw-sh)
          } else {
            ET_tot <- 0
          } 
          S_tot <- S_tot-ET_tot/(n*Zr)
          
          # Add Irrigation
          Q_int <- Q_day/interval
          if(S_tot<smin) {
            I_dm <- (smax-S_tot)*(n*Zr) # m
            I_tot <- min(I_dm,Q_int/A_tot)
            Q_int <- Q_int-A_tot*I_tot
          } else {
            I_tot <- 0
          } 
          S_tot <- S_tot+I_tot/(n*Zr) 
          
          # Save Interval Values
          S_tot.temp[k] <- S_tot
          I_tot.temp[k] <- I_tot
          ET_tot.temp[k] <- ET_tot
          R_tot.temp[k] <- R_tot 
          L_tot.temp[k] <- L_tot 
        }
        
        # Save Daily Values
        S_tot.v[((j-1)*interval+1):(j*interval)] <- S_tot.temp  
        I_tot.v[j] <- sum(I_tot.temp)
        ET_tot.v[j] <- sum(ET_tot.temp)
        R_tot.v[j] <- sum(R_tot.temp)
        L_tot.v[j] <- sum(L_tot.temp)
        c_day.v[j] <- c_day
        
        # Update Iterative Data Structures
        S_tot.old <- S_tot.v[(j*k)]
      }
      
      # Calculate Average Static Stress
      Z_tot.v <- (ss-S_tot.v[(interval*buffer+1):(interval*(Tseas+buffer))])/
        (ss-sw) # dim
      Z_tot.v <- Z_tot.v[Z_tot.v>0]
      Z_tot.v[Z_tot.v>1] <- 1
      if(length(Z_tot.v)>0) {
        Z_sub <- mean(Z_tot.v^q_sub) # dim
        Z_mar <- mean(Z_tot.v^q_mar) # dim
      } else {
        Z_sub <- 0 # dim
        Z_mar <- 0 # dim
      }
      
      # Calculate Crossing Parameters
      S_uns.v <- which(S_tot.v[(interval*buffer+1):(interval*(Tseas+buffer))]>=ss)
      Tss.v <- diff(c(0,S_uns.v,(interval*Tseas+1)))-1
      Tss.v <- Tss.v[Tss.v>0]
      nss <- length(Tss.v) # dim
      if(nss>0) {
        Tss <- mean(Tss.v)/interval # days
      } else {
        Tss <- 0
      }
      
      # Calculate Dynamic Stress
      theta_sub <- ((Z_sub*Tss)/(k_sub*Tseas))^(nss^-r_sub)
      if(theta_sub>1) {theta_sub <- 1} # dim
      theta_mar <- ((Z_mar*Tss)/(k_mar*Tseas))^(nss^-r_mar)
      if(theta_mar>1) {theta_mar <- 1} # dim
      
      # Calculate Crop Yields
      Y_sub <- A_sub*Ymax_sub*(1-theta_sub) # kg
      Y_mar <- A_mar*Ymax_mar*(1-theta_mar) # kg
      
      # Calculate total volume of water applied
      V_tot <- A_tot*sum(I_tot.v[(buffer+1):(Tseas+buffer)]) # m
      
      # Calculate Total Subsistence Crop Available
      Y_tot <- Y_sub+(cy_mar*Y_mar-sum(c_day.v[(buffer+1):(Tseas+buffer)]))/
        cy_sub # kg
      
      # Calculate Total Profit
      G_tot <- cy_sub*Y_sub+cy_mar*Y_mar-cw # $
      if(G_tot<=0) {G_tot <- 0}
      
      # Update Global Data Structures
      S_tot.m[,i] <- S_tot.v[(interval*buffer+1):(interval*(Tseas+buffer))]  
      I_tot.m[,i] <- I_tot.v[(buffer+1):(Tseas+buffer)]
      ET_tot.m[,i] <- ET_tot.v[(buffer+1):(Tseas+buffer)]
      R_tot.m[,i] <- R_tot.v[(buffer+1):(Tseas+buffer)]
      L_tot.m[,i] <- L_tot.v[(buffer+1):(Tseas+buffer)]
      theta_sub.v[i] <- theta_sub
      theta_mar.v[i] <- theta_mar
      Y_sub.v[i] <- Y_sub
      Y_mar.v[i] <- Y_mar
      Y_tot.v[i] <- Y_tot
      G_tot.v[i] <- G_tot
    }
    
    ##Update Superglobal Data Structures
    S_tot.a2[,,h] <- S_tot.m 
    S_tot_dens.l2[[h]] <- density(c(S_tot.m))
    S_tot_dens.m[1,h] <- min(S_tot_dens.l2[[h]]$x)
    S_tot_dens.m[2,h] <- max(S_tot_dens.l2[[h]]$x)
    S_tot_dens.m[3,h] <- max(S_tot_dens.l2[[h]]$y)
    I_tot.a2[,,h] <- I_tot.m
    I_tot_dens.l2[[h]] <- density(I_tot.m)
    I_tot_dens.m[1,h] <- min(I_tot_dens.l2[[h]]$x)
    I_tot_dens.m[2,h] <- max(I_tot_dens.l2[[h]]$x)
    I_tot_dens.m[3,h] <- max(I_tot_dens.l2[[h]]$y)
    Y_sub.m[,h] <- Y_sub.v
    Y_mar.m[,h] <- Y_mar.v
    Y_sub_dens.l2[[h]] <- density(Y_sub.v)
    Y_mar_dens.l2[[h]] <- density(Y_mar.v)
    Y_com_dens.m[1,h] <- min(Y_sub_dens.l2[[h]]$x,Y_mar_dens.l2[[h]]$x)
    Y_com_dens.m[2,h] <- max(Y_sub_dens.l2[[h]]$x,Y_mar_dens.l2[[h]]$x)
    Y_com_dens.m[3,h] <- max(Y_sub_dens.l2[[h]]$y,Y_mar_dens.l2[[h]]$y)
    Y_tot.m[,h] <- Y_tot.v
    Y_tot_dens.l2[[h]] <- density(Y_tot.v)
    Y_tot_dens.m[1,h] <- min(Y_tot_dens.l2[[h]]$x)
    Y_tot_dens.m[2,h] <- max(Y_tot_dens.l2[[h]]$x)
    Y_tot_dens.m[3,h] <- max(Y_tot_dens.l2[[h]]$y)
    G_tot.m[,h] <- G_tot.v
    G_tot_dens.l2[[h]] <- density(G_tot.v,na.rm=TRUE)
    G_tot_dens.m[1,h] <- min(G_tot_dens.l2[[h]]$x)
    G_tot_dens.m[2,h] <- max(G_tot_dens.l2[[h]]$x)
    G_tot_dens.m[3,h] <- max(G_tot_dens.l2[[h]]$y)
  }
  
  # Update Superduperglobal Data Structures
  S_tot.a1[,,,g] <- S_tot.a2
  S_tot_dens.l1[[g]] <- S_tot_dens.l2
  S_tot_dens.a[,,g] <- S_tot_dens.m
  I_tot.a1[,,,g] <- I_tot.a2
  I_tot_dens.l1[[g]] <- I_tot_dens.l2
  I_tot_dens.a[,,g] <- S_tot_dens.m
  Y_sub.a[,,g] <- Y_sub.m
  Y_mar.a[,,g] <- Y_mar.m
  Y_sub_dens.l1[[g]] <- Y_sub_dens.l2  
  Y_mar_dens.l1[[g]] <- Y_mar_dens.l2 
  Y_com_dens.a[,,g] <- Y_com_dens.m
  Y_tot.a[,,g] <- Y_tot.m
  Y_tot_dens.l1[[g]] <- Y_tot_dens.l2  
  Y_tot_dens.a[,,g] <- Y_tot_dens.m
  G_tot.a[,,g] <- G_tot.m
  G_tot_dens.l1[[g]] <- G_tot_dens.l2  
  G_tot_dens.a[,,g] <- G_tot_dens.m
}

# Find Rainfall Percentiles
P_tot.m <- P_tot.m[(interval*buffer+1):(interval*(Tseas+buffer)),]
P_sea.v <- colSums(P_tot.m)
perc <- ecdf(P_sea.v)
P_perc.v <- perc(P_sea.v)
P_10.v <- which(P_perc.v>=0.05 & P_perc.v<=0.15)
P_50.v <- which(P_perc.v>=0.45 & P_perc.v<=0.55)
P_90.v <- which(P_perc.v>=0.85 & P_perc.v<=0.95)

# Find Yield and Return Quantiles
stat <- function(x) {
  return(length(x[x>Ymin])/length(x))
}

stat_tot.a[1,,] <- apply(Y_tot.a,c(2,3),function(x){return(length(x[x>Ymin])/length(x))})
stat_tot.a[2,,] <- apply(G_tot.a,c(2,3),mean,na.rm=TRUE)
Y_10.a <- Y_tot.a[P_10.v,,]
G_10.a <- G_tot.a[P_10.v,,]
stat_10.a[1,,] <- apply(Y_10.a,c(2,3),function(x){return(length(x[x>Ymin])/length(x))})
stat_10.a[2,,] <- apply(G_10.a,c(2,3),mean,na.rm=TRUE)
Y_50.a <- Y_tot.a[P_50.v,,]
G_50.a <- G_tot.a[P_50.v,,]
stat_50.a[1,,] <- apply(Y_50.a,c(2,3),function(x){return(length(x[x>Ymin])/length(x))})
stat_50.a[2,,] <- apply(G_50.a,c(2,3),mean,na.rm=TRUE)
Y_90.a <- Y_tot.a[P_90.v,,]
G_90.a <- G_tot.a[P_90.v,,]
stat_90.a[1,,] <- apply(Y_90.a,c(2,3),function(x){return(length(x[x>Ymin])/length(x))})
stat_90.a[2,,] <- apply(G_90.a,c(2,3),mean,na.rm=TRUE)

stat.a2 <- abind(stat_tot.a,stat_10.a,stat_50.a,stat_90.a,along=4)

# Define risk aversion
risk.v <- Ymin*cy_sub*c(0, 0.4, 0.8, 1.2)
stat2.a2 <- array(NA, dim = c(length(risk.v), length(frac_sub.v), 
                              length(var.v), 4))
for(i in 1:length(risk.v))  {
  stat2.a2[i , , , ] <- stat.a2[2 , , , ] + risk.v[i] * stat.a2[1 , , , ]
}


## Plots!
# Set Working Directory
setwd(paste0(getwd(), "/Results"))
areas <- paste0(as.character(frac_sub.v * 100), "%")
irr_steps <- c("No Irrigation", "Irrigation Every Third Day", 
               "Irrigation Every Day")
mag_steps <- paste0("Qrte =  ", as.character(var.v), " m^3/day") 
cw_steps <- paste0("Cw =  $", as.character(var.v))
rat_steps <- paste0("Ratio = ", as.character(var.v))
steps <- irr_steps
patterns <- c(1, 2, 4, 5, 6)
partitions <- c("All Seasons", "Seasons in the 10th Quantile of Rainfall", 
                "Seasons in the 50th Quantile of Rainfall", 
                "Seasons in the 90th Quantile of Rainfall")
values <- paste0("Value of Successful Season =  $", as.character(risk.v))

rec <- sample(seq(1, trials), 1)
jpeg(filename = "smtrace.jpg", height = 500, width = 1000, quality = 100,
     bg = "white")
#quartz(width = 8.5, height = 5)
par(mar = c(4.5, 4.5, 1, 1) + 0.1)
plot(S_tot.a1[, rec, 3, 1], type = 'l', lwd = 2, col = "black", lty = 1, 
     cex.lab = 1.6, cex.axis = 1.5, font.main = 1, xlim = c(0, interval * Tseas),
     ylim = c(0, 1), xaxt = "n", xlab = "Days", ylab = "Soil Moisture (dim)", 
     main = NA)
for(i in 2:length(var.v)){
  lines(S_tot.a1[, rec, 3, i], type = 'l', lwd = 2, lty = patterns[i])
}
abline(h = sw, lty = 3)
abline(h = ss, lty = 3)
abline(h = sfc, lty = 3)
axis(1, at = seq(0, (interval * Tseas), (interval * Tseas / 5)), 
     labels = paste(seq(0, Tseas, (Tseas / 5))), cex.axis = 1.5)
legend("topright", legend = steps, lwd = 2, lty = patterns[1:length(var.v)], 
       cex = 1.5, bg = "white")
dev.off()

jpeg(filename="yperc.jpg", height = 700, width = 1250, quality = 100, 
     bg = "white")
#quartz(width = 8.5, height = 5)
layout(mat = matrix(1:4, nrow = 2, byrow = TRUE), height = 1)
par(mar = c(4.5, 4.5, 3, 1) + 0.1)
for(i in 1:4) {
  plot(100 * frac_sub.v, stat.a2[1, , 1, i], type = 'b', col = "black", lwd = 2, 
       cex.lab = 1.75, cex.axis = 1.5, cex.main = 2, ylim = c(0, 1), 
       xlab = "Percent Area Planted with Subsistence Crop", 
       ylab = "Subsistence Threshold Success Rate", main = partitions[i])
  for(j in 2:length(var.v)){
    lines(100 * frac_sub.v, stat.a2[1, , j, i], type = 'b', col = "black", 
          lwd = 2, lty = patterns[j])
  }
  legend("bottomright", legend = steps, lwd = 2, lty = patterns[1:length(var.v)], 
         cex = 1.75, bg = "white")
}
dev.off()

jpeg(filename = "gav.jpg", height = 700, width = 1250, quality = 100, 
     bg = "white")
#quartz(width = 8.5, height = 5)
layout(mat = matrix(1:4, nrow = 2, byrow = TRUE), height = 1)
par(mar=c(4.5,4.5,3,1)+0.1)
for(i in 1:4) {
  plot(100 * frac_sub.v, stat.a2[2, , 1, i], type = 'b', col = "black", lwd = 2, 
       lty = patterns[1], cex.lab = 1.75, cex.axis = 1.5, cex.main = 2, 
       ylim = c(0, max(1.1 * stat.a2[2, , , ], na.rm = TRUE)), 
       xlab = "Percent Area Planted with Subsistence Crop", 
       ylab = "Average Net Return ($)", main = partitions[i]) 
  for(j in 2:length(var.v)){
    lines(100 * frac_sub.v, stat.a2[2, , j, i], type = 'b', col = "black", 
          lwd = 2, lty = patterns[j])
  }
  legend("bottomright", legend = steps, lwd = 2, lty = patterns[1:length(var.v)], 
         cex = 1.75, bg = "white")
}
dev.off()

jpeg(filename = "tot.jpg", height = 700, width = 1250, quality = 100, 
     bg = "white")
# quartz(width = 8.5, height = 5)
layout(mat = matrix(1:4, nrow = 2, byrow = TRUE), height = 1)
par(mar = c(4.5, 4.5, 3, 1) + 0.1)
for(i in 1:4) {
  plot(100 * frac_sub.v, stat2.a2[i, , 2, 3], type = 'b', col = "black", lwd = 2, 
       lty=patterns[1], cex.lab = 1.75, cex.axis = 1.5, cex.main = 2, 
       xlab = "Percent Area Planted with Subsistence Crop", 
       ylab = "Average Total Value ($)", main = values[i]) 
}
dev.off()


## End time counter
elptm <- proc.time() - ptm