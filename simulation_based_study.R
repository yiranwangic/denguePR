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
library(dplyr)

# calculate the average temperature for the day t of one year
df <- read.csv('San_Juan_Temperature_Data.csv')
df$date <- seq(as.Date("1956-1-1"),by="day",length.out=nrow(df))
df <- df[df$date > "1990-4-29", ]
df <- df[df$date < "2013-4-30", ]
df$date_1 <- format(df$date, format="%m-%d")
doy_temperature_average <- ddply(df, ~date_1, colwise(mean, ~TAVG))
# remove the 29th Feb in leap year
doy_temperature_average <- doy_temperature_average[-60,] 
# simulate 1002-year temperature data
doy_temperature_average<-as.data.frame(doy_temperature_average)
temp <- rep(doy_temperature_average$TAVG,1002)

# temperature-dependent bt
a <- 0.0015  # randomly assigned value governing the height of thermal response curve 
Tmin <- 17.8 # values from mordecai's paper governing the shape of the thermal response curve
Tmax <- 34.6
Topt <- 29.1
m <- Topt*(Tmin-Topt)/(2*(Topt^2)-2*Topt*Tmax-Tmin*Topt+Tmin*Tmax)
B <- a*temp*((temp-Tmin))*(Tmax-temp)^(1/m)

#--- Create compartments ---#
nT <- 365*1002-1 # N time points (days)
nt2 <- nT/7
indl <- seq(1,nT,7)
indu <- seq(7,nT+6,7)
ind <- sort(rep(seq(1,nt2),7))
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
kappa <- 1/6 # human-mosquito-human transmission time
sigma <- 1/4 # infectious period
pop <- 5000000 # population
gamma <- 1/80/365 # birth/death rate (assumed to be equal)
rho <- 0.1 # assumed reporting rate
delta <- 0.5 # reporting rate of primary infections (relative to rho)
phi <- 2 # ADE factor
zeta <- 1/0.7/365 # cross-protection period

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
  
  # weekly cases in the week that the day belongs to
  for(s in 1:4) cases[t,s] <- round(rho*omega*(sum(e2[indl[ind[t]]:indu[ind[t]],s]) + delta*sum(e1[indl[ind[t]]:indu[ind[t]],s]))*pop)
  casetot[t] = sum(cases[t,])
}


# check that population always sums to 1 over time
tot <- su + rowSums(i1)+ rowSums(i2)+ rowSums(e1)+ rowSums(e2)+ rowSums(moP)+rowSums(moS)+ r
plot(tot, type='l') 

# set state variables at t=361816 as the initial conditions for the second simulation
isu <- su[361816] 
ir <- r[361816] 
iv1 <- sum(e1[361816,]) + sum(i1[361816,])
imop <- moP[361816,]  
imos <- moS[361816,]
iv2 <- sum(e2[361816,]) + sum(i2[361816,])
pSERO <- cases[361816,]/casetot[361816]

# start the model again and run the simulation for 10 years
temp <- temp[361816:365462]
# temperature-dependent bt
a <- 0.0015  # randomly assigned value governing the height of thermal response curve 
Tmin <- 17.8 # values from mordecai's paper governing the shape of the thermal response curve
Tmax <- 34.6
Topt <- 29.1
m <- Topt*(Tmin-Topt)/(2*(Topt^2)-2*Topt*Tmax-Tmin*Topt+Tmin*Tmax)
B <- a*temp*((temp-Tmin))*(Tmax-temp)^(1/m)

#--- Create compartments ---#
nT <- length(temp) # N time points (days)
nt2 <- nT/7
indl <- seq(1,nT,7)
indu <- seq(7,nT+6,7)
ind <- sort(rep(seq(1,nt2),7))
su <- rep(0, nT) # susceptible
e1 <- matrix(0, nrow=nT, ncol=4) # exposed, 1st infection
e2 <- matrix(0, nrow=nT, ncol=4) # exposed, 2nd infection
i1 <- matrix(0, nrow=nT, ncol=4) # infectious, 1st infection
i2 <- matrix(0, nrow=nT, ncol=4) # infectious, 2nd infection
moP <- matrix(0, nrow=nT, ncol=4) # cross protected
moS <- matrix(0, nrow=nT, ncol=4) # monotypic
r <- rep(0, nT) # recovered/immune

v1 <- rep(0, nT) 
v2 <- rep(0, nT) 

lam <- matrix(0, nrow=nT, ncol=4) # transmission intensity
mz <- matrix(0, nrow=nT, ncol=4) # infectious mosquitoes
cases <- matrix(0, nrow=nt2, ncol=4) # reported cases
casetot <- rep(0,nt2)

#--- Parameters ---#
omega <- 1/6 # IIP
kappa <- 1/6 # human-mosquito-human transmission time
sigma <- 1/4 # infectious period
pop <- 5000000 # population
gamma <- 1/80/365 # birth/death rate (assumed to be equal)
rho <- 0.1 # assumed reporting rate
delta <- 0.5 # reporting rate of primary infections (relative to rho)
phi <- 2 # ADE factor
zeta <- 1/0.7/365 # cross-protection period

#--- Initial conditions ---#
su[1] <- isu # initial susceptible fraction
r[1] <- ir # initial recovered fraction
v1[1] <- iv1
moP[1,] <- imop 
moS[1,] <- imos
v2[1] <- iv2
tot <- su[1] + r[1] + v1[1] + sum(moP[1,]) + sum(moS[1,]) + v2[1]
tot

# approximated initial values for E and I generated by using the quasi-equillibrium assumption and pSERO
e1[1,] <- v1[1]*pSERO*(sigma/(sigma+omega))
i1[1,] <- v1[1]*pSERO*(omega/(sigma+omega))
e2[1,] <- v2[1]*pSERO*(sigma/(sigma+omega))
i2[1,] <- v2[1]*pSERO*(omega/(sigma+omega))

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
  
  # aggregate to weekly cases
  for(t in 1:nt2) for(s in 1:4) cases[t,s] <- round(rho*omega*(sum(e2[indl[t]:indu[t],s]) + delta*sum(e1[indl[t]:indu[t],s]))*pop)
  for(t in 1:nt2) casetot[t] = sum(cases[t,])
}

# plot simulated 10-year case data
matplot(su[t=(1:nT)], type='l', ylim=c(0,1), xlab='Time (Years)', ylab='Susceptible',xaxt='n')
matplot(r[t=(1:nT)], type='l', ylim=c(0,1), xlab='Time (Years)', ylab='Recovered',xaxt='n')


col_set <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
matplot(cases,type='l',xlab='Time (Years)',ylab='Serotype-specific cases',lty=1,col=col_set,xaxt='n',main='')
matplot(casetot, type='l', xlab='Time (Years)', ylab='Total cases',xaxt='n')

# input data
data <- list()
# remove the first two weeks where the equillibrum was skewed by the second initiation of simulation
data$Nt <- nT-14
data$Nt2 <- data$Nt/7
data$indL <- seq(1,data$Nt,7)
data$indU <- seq(7,data$Nt+6,7)

data$Ncases <- cases[-(1:2),]
data$NcasesTot <- casetot[-(1:2)]

# input parameters
data$pop <- pop
data$gamma <- gamma
data$sigma <- sigma
data$omega <- omega
data$kappa <- kappa
data$temp <- temp[-(1:14)]
data$pSERO <- pSERO
data$rho <- rho
data$delta <- delta
data$phi <- phi
data$zeta <- zeta
data$a <- a
data$Tmin <- Tmin
data$Tmax <- Tmax
data$m <- m
data$alpha <- c(20,60,0.4,rep(5,8),0.4)/5

# the actual values of initial conditions
su[14]
r[14]
sum(e1[14,]) + sum(i1[14,])
moP[14,]
moS[14]
sum(e2[14,]) + sum(i2[14,])

# fit the model to simulated data 
alpha <- c(20,60,0.4,rep(5,8),0.4)/5
initI <- rdirichlet(100,alpha)
ii1 <- as.vector(initI[1,])
ii2 <- as.vector(initI[2,])
ii3 <- as.vector(initI[3,])
inits_chain1 <- list(init=ii1, log_a=-6.5, Tmin=17, Tmax=35)  
inits_chain2 <- list(init=ii2, log_a=-6.4, Tmin=16, Tmax=34)
inits_chain3 <- list(init=ii3, log_a=-6.6, Tmin=18, Tmax=36)
check_cmdstan_toolchain(fix=T)
set_cmdstan_path()
mod <- cmdstan_model('simulation_based_study.stan', pedantic=T)
fit <- mod$sample(data=data, chains=3, parallel_chains=3, iter_sampling=1000, refresh=10, iter_warmup=2000, save_warmup=FALSE, init=list(inits_chain1,inits_chain2,inits_chain3))
stanfit <- rstan::read_stan_csv(fit$output_files())
chains <- rstan::extract(stanfit)

# Check convergence
traceplot1 <- traceplot(stanfit, pars=c('init','lp__'))
traceplot1
traceplot2 <- traceplot(stanfit, pars=c('a','Tmin','Tmax'))
traceplot2

# check model fitting
fitS <- data.frame(t=rep(seq(1:data$Nt2),4),sero=sort(rep(seq(1,4),data$Nt2)),obs=NA, med=NA,ciL=NA,ciU=NA)
fitTot <- data.frame(t=seq(1,data$Nt2),obs=data$NcasesTot, med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt2){
  for(s in 1:4){
    fitS[fitS$t==i & fitS$sero==s,4:6] <- quantile(chains$predCases[,i,s], c(0.5,0.025,0.975))
    fitS[fitS$t==i & fitS$sero==s,3] <- data$Ncases[i,s]
  }
}
for(i in 1:data$Nt2){
  fitTot[fitTot$t==i,3:5] <- quantile(chains$predCasesTot[,i], c(0.5,0.025,0.975))
  fitTot[fitTot$t==i,2] <- data$NcasesTot[i]
}
fitPlotS <- ggplot(fitS, aes(t,obs))+ geom_point(alpha=0.4, col='grey40')+
  theme_minimal()+ theme(legend.position='none')+ theme(text=element_text(size=16))+
  geom_line(aes(t,med,col=factor(sero)))+ geom_ribbon(aes(t, ymin=ciL,ymax=ciU,fill=factor(sero)),alpha=0.3)+
  xlab('Time')+ ylab('Cases')+ theme(legend.position='none')+ facet_wrap(~sero,scale='free_y')
fitPlotS
fitPlotTot <- ggplot(fitTot, aes(t,obs))+ geom_point(alpha=0.4, col='grey40')+
  theme_minimal()+ geom_line(aes(t,med), col='blue')+ xlab('Time')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='dodgerblue', alpha=0.5)+ ylab('Cases')+
  theme(text=element_text(size=16))
fitPlotTot

# plot model estimates vs actual values
fit1 <- as.data.frame(summary(stanfit,pars=c('a','Tmin','Tmax'), probs = c(0.025, 0.975))$summary)
fit1[,"real_value"] <- c(a,Tmin,Tmax)
fit1 <- data.frame(parameter=c('a','Tmin','Tmax'), real_value=fit1$real_value, mean=fit1$mean, lci=fit1$`2.5%`,uci=fit1$`97.5%`)

fit1_1 <- fit1[1,]
plotfit1_1 <- ggplot(fit1_1)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_1
fit1_2 <- fit1[2,]
plotfit1_2 <- ggplot(fit1_2)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_2
fit1_3 <- fit1[3,]
plotfit1_3 <- ggplot(fit1_3)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_3

fit2 <- as.data.frame(summary(stanfit,pars=c('init'), probs = c(0.025, 0.975))$summary)
fit2[,"real_value"] <- c(su[14], r[14], 
                         sum(e1[14,])+sum(i1[14,]),
                         moP[14,1], moP[14,2], moP[14,3], moP[14,4], 
                         moS[14,1], moS[14,2], moS[14,3], moS[14,4],
                         sum(e2[14,])+sum(i2[14,]))
fit2 <- data.frame(initial_state=c('S','R','V1','moP1','moP2','moP3','moP4','moS1','moS2','moS3','moS4','V2'), 
                   real_value=fit2$real_value, mean=(fit2$mean), 
                   lci=fit2$`2.5%`,uci=fit2$`97.5%`)

fit2_1 <- fit2[1,]
plotfit2_1 <- ggplot(fit2_1)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_1
fit2_2 <- fit2[2,]
plotfit2_2 <- ggplot(fit2_2)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_2
fit2_3 <- fit2[3,]
plotfit2_3 <- ggplot(fit2_3)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_3
fit2_4 <- fit2[4:7,]
plotfit2_4 <- ggplot(fit2_4)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_4
fit2_5 <- fit2[8:11,]
plotfit2_5 <- ggplot(fit2_5)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_5
fit2_6 <- fit2[12,]
plotfit2_6 <- ggplot(fit2_6)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_6

fit2 <- data.frame(summary(stanfit,pars=c('init'), probs = c(0.025, 0.975))$summary)
fit2[13,] <- data$pSERO[1]*(fit2[3,]+fit2[12,])+fit2[4,]+fit2[8,]
fit2[14,] <- data$pSERO[2]*(fit2[3,]+fit2[12,])+fit2[5,]+fit2[9,]
fit2[15,] <- data$pSERO[3]*(fit2[3,]+fit2[12,])+fit2[6,]+fit2[10,]
fit2[16,] <- data$pSERO[4]*(fit2[3,]+fit2[12,])+fit2[7,]+fit2[11,]
fit3 <- fit2[13:16,]
fit3[5,] <- fit2[1,]
fit3[6,] <- fit2[2,]
fit3[,"real_value"] <- c((e1[1,1]+i1[1,1]+moP[1,1]+moS[1,1]+e2[1,1]+i2[1,1]), 
                         (e1[1,2]+i1[1,2]+moP[1,2]+moS[1,2]+e2[1,2]+i2[1,2]),
                         (e1[1,3]+i1[1,3]+moP[1,3]+moS[1,3]+e2[1,3]+i2[1,3]), 
                         (e1[1,4]+i1[1,4]+moP[1,4]+moS[1,4]+e2[1,4]+i2[1,4]),
                         su[1], r[1])
fit3[,"initial_state"] <- c(c('I_1','I_2','I_3','I_4','S','R'))
fit3 <- data.frame(initial_state=c('I_1','I_2','I_3','I_4','S','R'), 
                   real_value=fit3$real_value, mean=fit3$mean, 
                   lci=fit3$`X2.5.`,uci=fit3$`X97.5.`)
fit3_1 <- fit3[1:4,]
plotfit3_1 <- ggplot(fit3_1)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit3_1


# save results
png(filename='traceplot1.png', width=20, height=14, res=400, units='cm')
plot(traceplot2)
dev.off()
png(filename='traceplot2.png', width=20, height=14, res=400, units='cm')
plot(traceplot3)
dev.off()
png(filename='FitTotalCases.png', width=20, height=14, res=400, units='cm')
plot(fitPlotTot)
dev.off()
png(filename='FitSeroCases.png', width=20, height=14, res=400, units='cm')
plot(fitPlotS)
dev.off()
png(filename='1.png', width=20, height=20, res=400, units='cm')
plot(plotfit1_1)
dev.off()
png(filename='2.png', width=20, height=20, res=400, units='cm')
plot(plotfit1_2)
dev.off()
png(filename='3.png', width=20, height=20, res=400, units='cm')
plot(plotfit1_3)
dev.off()
png(filename='4.png', width=20, height=20, res=400, units='cm')
plot(plotfit1_4)
dev.off()
png(filename='5.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_1)
dev.off()
png(filename='6.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_2)
dev.off()
png(filename='7.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_3)
dev.off()
png(filename='8.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_4)
dev.off()
png(filename='9.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_5)
dev.off()
png(filename='10.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_6)
dev.off()
png(filename='11.png', width=20, height=20, res=400, units='cm')
plot(plotfit3_1)
dev.off()
