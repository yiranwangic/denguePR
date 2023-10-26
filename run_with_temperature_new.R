#----- Simulate dengue serotype model -----#
# assumes complete immunity after 2 infections
# varying serotype proportions
setwd('/Users/wangyr/Desktop')
library(ggplot2)
library(plyr)
library(cmdstanr)
library(rstan)
library(data.table)
library(ggplot2)
library(DirichletReg)
library(posterior)


# calculate the average temperature for the day t of one year
df <- read.csv('San_Juan_Temperature_Data.csv')
df$date <- seq(as.Date("1956-1-1"),by="day",length.out=nrow(df))
df <- df[df$date > "1990-4-29", ]
df <- df[df$date < "2013-4-30", ]
df$date_1 <- format(df$date, format="%m-%d")
doy_temperature_average <- ddply(df, ~date_1, colwise(mean, ~TAVG))
# remove the 29th Feb in leap year
doy_temperature_average <- doy_temperature_average[-60,] 
# simulate 1010-year temperature data
doy_temperature_average<-as.data.frame(doy_temperature_average)
temp <- rep(doy_temperature_average$TAVG,1010)

# temperature-dependent beta_t - general briere function
a <- 0.0015
Tmin <- 17.8 # values from Mordecai's paper
Tmax <- 34.6
Topt <- 29.1
m <- Topt*(Tmin-Topt)/(2*(Topt^2)-2*Topt*Tmax-Tmin*Topt+Tmin*Tmax)

B <- a*temp*((temp-Tmin))*((Tmax-temp)^(1/m))
B <- ifelse(temp<Tmin,0,ifelse(temp>Tmax,0,B))


#--- Create compartments ---#

nT <- 365*1010 # N time points (days)
su <- rep(0, nT) # susceptible
e1 <- matrix(0, nrow=nT, ncol=4) # exposed, 1st infection
e2 <- matrix(0, nrow=nT, ncol=4) # exposed, 2nd infection
i1 <- matrix(0, nrow=nT, ncol=4) # infectious, 1st infection
i2 <- matrix(0, nrow=nT, ncol=4) # infectious, 2nd infection
moP <- matrix(0, nrow=nT, ncol=4) # cross protected
moS <- matrix(0, nrow=nT, ncol=4) # monotypic
r <- rep(0, nT) # recovered/immune
lam <- matrix(0, nrow=nT, ncol=4) # transmission intensity
mz <- matrix(0, nrow=nT, ncol=4) # infectious mosquitoes
cases <- matrix(0, nrow=nT, ncol=4) # reported cases
casetot <- rep(0,nT)

#--- Parameters ---#
omega <- 1/6 # IIP
kappa <- 1/6 # EIP
sigma <- 1/4 # infectious period
pop <- 10000000 # population of 10 million
gamma <- 1/80/365 # birth/death rate (assumed to be equal)
rho <- 0.1 # assumed reporting rate
delta <- 0.05 # reporting rate of primary infections (relative to rho)
phi <- 2 # relative increased transmissability of 2nd infections
zeta <- 1/0.7/365 # TCI period

#--- Initial conditions ---#
inf <- c(0.000003,0.000005,0.000001,0.000008)
for(s in 1:4) i1[1,s] <- inf[s] 
for(s in 1:4) i2[1,s] <- inf[s]
for(s in 1:4) moP[1,s] <- 0  
for(s in 1:4) moS[1,s] <- 0
r[1] <- 0 # initial recovered fraction
su[1] <- 1 - (sum(i1[1,])+ sum(i2[1,])+ sum(moP[1,])+ sum(moS[1,])+ sum(r[1])) # initial susceptible fraction
mz[1,s] = B[1]*(i1[1,s]+phi*i2[1,s]) # initial infectious mosquito fraction

#--- Simulate ---#
for(t in 1:(nT-1)){
  
  for(s in 1:4) mz[t+1,s] <- mz[t,s]+ B[t]*sum(i1[t,s] + phi*i2[t,s]) - kappa*mz[t,s]
  for(s in 1:4) lam[t,s] <- kappa*mz[t,s]
  su[t+1] <- su[t] - su[t]*sum(lam[t,]) - gamma*su[t] + gamma
  for(s in 1:4) e1[t+1,s] <- e1[t,s] + su[t]*lam[t,s] - omega*e1[t,s] - gamma*e1[t,s]
  for(s in 1:4) i1[t+1,s] <- i1[t,s] + omega*e1[t,s] - sigma*i1[t,s] - gamma*i1[t,s]
  for(s in 1:4) moP[t+1,s] <- moP[t,s] + sigma*i1[t,s] - moP[t,s]*zeta - gamma*moP[t,s]
  for(s in 1:4) moS[t+1,s] <- moS[t,s] + moP[t,s]*zeta - moS[t,s]*sum(lam[t,-s]) - gamma*moS[t,s]
  for(s in 1:4) e2[t+1,s] <- e2[t,s] + sum(moS[t,-s])*lam[t,s] - omega*e2[t,s] - gamma*e2[t,s]
  for(s in 1:4) i2[t+1,s] <- i2[t,s] + omega*e2[t,s] - sigma*i2[t,s] - gamma*i2[t,s]
  r[t+1] <- r[t] + sigma*sum(i2[t,]) - gamma*r[t]
  
  for(s in 1:4) cases[t,s] <- rho*omega*(e2[t,s] + delta*e1[t,s])*pop
  casetot[t] = sum(cases[t,])
  if(t %in% seq(1:1010)*365) print(paste(t))
}


# check that population always sums to 1 over time
tot <- su + rowSums(i1)+ rowSums(i2)+ rowSums(e1)+ rowSums(e2)+ rowSums(moP)+rowSums(moS)+ r
plot(tot, type='l') 

# plot simulated case data
col_set <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
matplot(cases[t=(361816:368650),],type='l',xlab='Time (Years)',ylab='Serotype-specific cases',lty=1,col=col_set,xaxt='n',main='')
axis(1,at=((0:20)*5*365),labels=((0:20)*5),cex.axis=0.8)
matplot(casetot[t=(361350:368650)], type='l', xlab='Time (Years)', ylab='Total cases',xaxt='n')
axis(1,at=((0:20)*5*365),labels=((0:20)*5),cex.axis=0.8)


# create 20-year simulated data (choosing a starting point that with cases for all four serotypes)
data <- list()
data$Ncases <- round(cases[t=(361816:368650),])
data$NcasesTot <- rowSums(data$Ncases)
data$Nt <- nrow(data$Ncases)
data$pop <- 10000000
data$gamma <- 1/80/365
data$sigma <- 1/4
data$omega <- 1/6
data$kappa <- 1/6
data$alpha <- c(20,60,0.4,rep(5,4),rep(5,4),0.4)/2
data$temp <- temp[361816:368650]
data$pSERO <- round(cases[361816,])/sum(round(cases[361816,]))
data$delta <- delta
data$phi <- phi
data$zeta <- zeta
data$rho <- rho
data$m <- m
data$a <- a
data$Tmin <- Tmin
data$Tmax <- Tmax

# fit to the simulated data ~ poisson
alpha <- c(20,60,0.4,rep(5,4),rep(5,4),0.4)/1
initI <- rdirichlet(100,alpha)
ii1 <- as.vector(initI[1,])
ii2 <- as.vector(initI[2,])
ii3 <- as.vector(initI[3,])
inits_chain1 <- list(init=ii1)  
inits_chain2 <- list(init=ii2)
inits_chain3 <- list(init=ii3)
check_cmdstan_toolchain(fix=T)
set_cmdstan_path()
mod <- cmdstan_model('with_temperature_new.stan', pedantic=T)
fit <- mod$sample(data=data, chains=3, parallel_chains=3, iter_sampling=2000, refresh=10, iter_warmup=3000, save_warmup=TRUE, init=list(inits_chain1,inits_chain2,inits_chain3))
stanfit <- rstan::read_stan_csv(fit$output_files())
