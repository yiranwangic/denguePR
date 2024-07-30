#---- Fit SEIR model to San Juan serotype case data -----#
# libraries
library(cmdstanr)
library(rstan)
library(data.table)
library(ggplot2)
library(DirichletReg)

# read in Puerto Rico data
setwd("/Volumes/yw121")
df <- read.csv('san_juan_dengue_data.csv')
dates <- unique(df$week_start_date)
Ncases <- matrix(NA, nrow=length(dates), ncol=4)
for(s in 1:4) Ncases[,s] <- df[,3+s]

## Read in temp data
dft <- read.csv('~/Desktop/dengue/temp.csv')
# discard the last one/two days of one year to accommodate to the case data
df2 <- subset(dft,MM == 12 & DD == 31)
dft <- dft[-c(246,611,976,977,1342,1707,2072,2437,2438,2803,3168,3533,3898,3899,4264,
              4629,4994,5359,5360,5725,6090,6455,6820,6821,7186,7551,7916,8281,8282),]

# inputs for model
data <- list()
data$Ncases <- Ncases[-(1:18),] # choose the starting point with cases for DENV1, DENV2 and DENV4
data$NcasesTot <- df$total_cases[-(1:18)]
data$Nt <- nrow(data$Ncases)
data$Nt2 <- data$Nt*7
# population for each year of the time-series
data$pop <- c(rep(2327000,17),
              rep(2345000,52),
              rep(2363000,52),
              rep(2381000,52),
              rep(2400000,52),
              rep(2418000,52),
              rep(2437000,52),
              rep(2456000,52),
              rep(2475000,52),
              rep(2494000,52),
              rep(2508000,52),
              rep(2505000,52),
              rep(2502000,52),
              rep(2499000,52),
              rep(2496000,52),
              rep(2493000,52),
              rep(2490000,52),
              rep(2487000,52),
              rep(2484000,52),
              rep(2481000,52),
              rep(2478000,52),
              rep(2475000,52),
              rep(2472000,52),
              rep(2469000,17))
data$sigma <- 1/4
data$omega <- 1/6
data$kappa <- 1/6
data$phi <- 2
data$zeta <- 1/0.7/365
data$m <- m
# annual average birth rate for each year of the time-series
data$br<- c(rep(4.96e-05,17*7),
            rep(4.90e-05,52*7),
            rep(4.84e-05,52*7),
            rep(4.78e-05,52*7),
            rep(4.67e-05,52*7),
            rep(4.55e-05,52*7),
            rep(4.44e-05,52*7),
            rep(4.32e-05,52*7),
            rep(4.21e-05,52*7),
            rep(4.13e-05,52*7),
            rep(4.06e-05,52*7),
            rep(3.98e-05,52*7),
            rep(3.91e-05,52*7),
            rep(3.84e-05,52*7),
            rep(3.77e-05,52*7),
            rep(3.71e-05,52*7),
            rep(3.64e-05,52*7),
            rep(3.57e-05,52*7),
            rep(3.51e-05,52*7),
            rep(3.39e-05,52*7),
            rep(3.27e-05,52*7),
            rep(3.15e-05,52*7),
            rep(3.04e-05,52*7),
            rep(2.92e-05,17*7))
# annual average exit rate for each year of the time-series
data$er<- c(rep(2.84e-05,17*7),
            rep(2.80e-05,52*7),
            rep(2.75e-05,52*7),
            rep(2.59e-05,52*7),
            rep(2.61e-05,52*7),
            rep(2.40e-05,52*7),
            rep(2.30e-05,52*7),
            rep(2.20e-05,52*7),
            rep(2.10e-05,52*7),
            rep(2.60e-05,52*7),
            rep(4.39e-05,52*7),
            rep(4.31e-05,52*7),
            rep(4.24e-05,52*7),
            rep(4.17e-05,52*7),
            rep(4.10e-05,52*7),
            rep(4.03e-05,52*7),
            rep(3.97e-05,52*7),
            rep(3.90e-05,52*7),
            rep(3.84e-05,52*7),
            rep(3.72e-05,52*7),
            rep(3.60e-05,52*7),
            rep(3.49e-05,52*7),
            rep(3.37e-05,52*7),
            rep(3.25e-05,17*7))
data$alpha <- c(300,4000,20,rep(30,4),rep(90,4))/100
data$alpha1 <- c(300,4000,rep(5,4),rep(30,4),rep(90,4))/100
indL <- seq(1,data$Nt*7,7)
indU <- seq(7,data$Nt*7+6,7)
data$indL <- indL
data$indU <- indU
data$ind <- sort(rep(seq(1,data$Nt),7))
data$pSample <- rowSums(data$Ncases)/data$NcasesTot
data$pSample[is.na(data$pSample)] <- 0
data$temp <- dft$TAVG[-(1:18*7)]
data$pSERO <- data$Ncases[1,]/sum(data$Ncases[1,])

# chain starting values
alpha <- c(300,4000,rep(5,4),rep(30,4),rep(90,4))/100
initI <- rdirichlet(100,alpha)
ii1 <- as.vector(initI[1,])
ii2 <- as.vector(initI[2,])
ii3 <- as.vector(initI[3,])
ii4 <- as.vector(initI[4,])
ii5 <- as.vector(initI[5,])
ii6 <- as.vector(initI[6,])
inits_chain1 <- list(init=ii1, log_a=-6.5, Tmin=17, Tmax=35, pv1=0.5, log_rho=-2, log_delta=-4.5, log_est=-14, phi=1.5, log_zeta=-6.1)  
inits_chain2 <- list(init=ii2, log_a=-6.4, Tmin=16, Tmax=34, pv1=0.4, log_rho=-2.1, log_delta=-4.4, log_est=-14.2, phi=1.6, log_zeta=-6.2)
inits_chain3 <- list(init=ii3, log_a=-6.6, Tmin=18, Tmax=36, pv1=0.6, log_rho=-2.2, log_delta=-4.3, log_est=-14.4, phi=1.7, log_zeta=-6.3)
inits_chain4 <- list(init=ii4, log_a=-6.3, Tmin=15, Tmax=33, pv1=0.3, log_rho=-2.3, log_delta=-4.2, log_est=-14.6, phi=1.8, log_zeta=-6.4)  
inits_chain5 <- list(init=ii5, log_a=-6.7, Tmin=19, Tmax=37, pv1=0.7, log_rho=-2.4, log_delta=-4.1, log_est=-13.8, phi=1.9, log_zeta=-6.5)
inits_chain6 <- list(init=ii6, log_a=-6.8, Tmin=20, Tmax=38, pv1=0.8, log_rho=-2.5, log_delta=-4, log_est=-13.6, phi=2, log_zeta=-6.6)
mod <- cmdstan_model('fit_to_san_juan_data14_2.stan', pedantic=T)

# fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path()
mod <- cmdstan_model('fit_to_san_juan_data14_2.stan', pedantic=T)
fit <- mod$sample(data=data, 
                  chains=6,
                  parallel_chain=6,
                  iter_sampling=1000,
                  refresh=10,
                  iter_warmup=2000, 
                  init=list(inits_chain1,inits_chain2,inits_chain3,inits_chain4, inits_chain5, inits_chain6))

