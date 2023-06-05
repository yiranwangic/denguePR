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
# simulate 1000-year temperature data
doy_temperature_average<-as.data.frame(doy_temperature_average)
temp <- rep(doy_temperature_average$TAVG,1000)

# temperature-dependent bt
a <- 0.000202*0.000491*4
Tmin <- 13.35
Tmax <- 37.46
B <- a*(temp^2)*((temp-Tmin)^2)*(Tmax-temp) 

#--- Create compartments ---#
nT <- 365000 # N time points (days)
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
pop <- 10000000 # population of 2.5 million
gamma <- 1/80/365 # birth/death rate (assumed to be equal)
rho <- 0.2 # assumed reporting rate
delta <- 0.05 # reporting rate of primary infections (relative to rho)
phi <- 1.5 # relative increased transmissability of 2nd infections
zeta <- 1/0.5/365 # TCI period

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
  if(t %in% seq(1:1000)*365) print(paste(t))
}


# check that population always sums to 1 over time
tot <- S + rowSums(I1)+ rowSums(I2)+ rowSums(E1)+ rowSums(E2)+ rowSums(MoP)+rowSums(MoS)+ R
plot(tot, type='l') 

# plot simulated case data
col_set <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
matplot(cases[t=(357701:365000),],type='l',xlab='Time (Years)',ylab='Serotype-specific cases',col=col_set,lty=1,xaxt='n')
axis(1,at=((0:20)*5*365),labels=((0:20)*5),cex.axis=0.8)
matplot(casetot[t=(357701:365000)], type='l', xlab='Time (Years)', ylab='Total cases',xaxt='n')
axis(1,at=((0:20)*5*365),labels=((0:20)*5),cex.axis=0.8)

# create the 20-year simulated data
data <- list()
data$Ncases <- round(cases[t=(357701:365000),])
data$NcasesTot <- round(casetot[t=(357701:365000)])
data$Nt <- 365*20
data$pop <- 10000000
data$gamma <- 1/80/365
data$sigma <- 1/4
data$omega <- 1/6
data$kappa <- 1/6
data$alpha <- c(20,60,rep(5.2,4),rep(5,4),rep(0.1,4))/2
data$temp <- rep(doy_temperature_average$TAVG,20)

# fit to the simulated data
ii <- c(0.21, 0.58-0.0002057, 1.6e-2, 3.44e-4, 1.37e-3, 2.5e-3, 3.75e-2,6.15e-2,5.11e-2,3.98e-2,1.24e-5,5.45e-5,1.68e-5,8e-6)
inits_chain1 <- list(log_rho=-1.88, init=ii, log_delta=-1.27,a=4.3e-7,Tmin=13.5, Tmax=37.3, phi=1.42,log_zeta=-5.2)  
inits_chain2 <- list(log_rho=-1.85, init=ii, log_delta=-1.3, a=4.1e-7, Tmin=13, Tmax=37.1, phi=1.5, log_zeta=-5.3)
inits_chain3 <- list(log_rho=-1.80, init=ii, log_delta=-1.25,a=4.2e-7,Tmin=13.3, Tmax=37.2, phi=1.3, log_zeta=-5.1)
check_cmdstan_toolchain(fix=T)
set_cmdstan_path()
mod <- cmdstan_model('simulation_based_study.stan', pedantic=T)
fit <- mod$sample(data=data, chains=3, parallel_chains=3, iter_sampling=800, refresh=10, iter_warmup=200, save_warmup=TRUE, init=list(inits_chain1,inits_chain2,inits_chain3))
stanfit <- rstan::read_stan_csv(fit$output_files())
chains <- rstan::extract(stanfit)

# Check convergence
traceplot1 <- traceplot(stanfit, pars=c('init','lp__'))
traceplot1
traceplot2 <- traceplot(stanfit, pars=c('rho','delta','a','Tmin','Tmax','phi','log_zeta'))
traceplot2

# Check parameter estimates vs actual values
fit1 <- as.data.frame(summary(stanfit,pars=c('rho','delta','a','Tmin','Tmax','phi','log_zeta'), probs = c(0.025, 0.975))$summary)
fit1[,"real_value"] <- c(rho,delta,a,Tmin,Tmax,phi,log(zeta)) 
fit1 <- data.frame(parameter=c('rho','delta','a','Tmin','Tmax','phi','log_zeta'), real_value=fit1$real_value, mean=fit1$mean, lci=fit1$`2.5%`,uci=fit1$`97.5%`)
fit1_1 <- fit1[1:2,]
plotfit1_1 <- ggplot(fit1_1)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_1
fit1_2 <- fit1[3,]
plotfit1_2 <- ggplot(fit1_2)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_2
fit1_3 <- fit1[4:5,]
plotfit1_3 <- ggplot(fit1_3)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_3
fit1_4 <- fit1[6,]
plotfit1_4 <- ggplot(fit1_4)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_4
fit1_5 <- fit1[7,]
plotfit1_5 <- ggplot(fit1_5)+
  geom_point(aes(x=parameter, y=mean),size=2)+ 
  geom_errorbar(aes(x=parameter, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=parameter,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit1_5

# Check initial condition estimates vs actual values
fit2 <- as.data.frame(summary(stanfit,pars=c('init'), probs = c(0.025, 0.975))$summary)
fit2[,"real_value"] <- c(su[357701], r[357701], (e1[357701,1]+i1[357701,1]+moP[357701,1]), (e1[357701,2]+i1[357701,2]+moP[357701,2]),
                         (e1[357701,3]+i1[357701,3]+moP[357701,3]), (e1[357701,4]+i1[357701,4]+moP[357701,4]), moS[357701,1], moS[357701,2],
                         moS[357701,3], moS[357701,4], (e2[357701,1]+i2[357701,1]), (e2[357701,2]+i2[357701,2]),
                         (e2[357701,3]+i2[357701,3]), (e2[357701,4]+i2[357701,4]))
fit2 <- data.frame(initial_state=c('S','R','V1_1','V1_2','V1_3','V1_4','moS1','moS2',
                                   'moS3','moS4','V2_1','V2_2','V2_3','V2_4'), 
                   real_value=fit2$real_value, mean=fit2$mean, 
                   lci=fit2$`2.5%`,uci=fit2$`97.5%`)

fit2_1 <- fit2[1:2,]
plotfit2_1 <- ggplot(fit2_1)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_1
fit2_2 <- fit2[3:6,]
plotfit2_2 <- ggplot(fit2_2)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_2
fit2_3 <- fit2[7:10,]
plotfit2_3 <- ggplot(fit2_3)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_3
fit2_4 <- fit2[11:14,]
plotfit2_4 <- ggplot(fit2_4)+
  geom_point(aes(x=initial_state, y=mean),size=2)+ 
  geom_errorbar(aes(x=initial_state, ymin=lci, ymax=uci, width=0.4))+theme(legend.position='none')+
  geom_point(aes(x=initial_state,y=real_value),col='red')+theme_classic()+theme(text=element_text(size=16))
plotfit2_4

# Time-series fit
fitS <- data.frame(t=rep(seq(1:data$Nt),4),sero=sort(rep(seq(1,4),data$Nt)),obs=NA, med=NA,ciL=NA,ciU=NA)
fitTot <- data.frame(t=seq(1,data$Nt),obs=data$NcasesTot, med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt){
  for(s in 1:4){
    fitS[fitS$t==i & fitS$sero==s,4:6] <- quantile(chains$predCases[,i,s], c(0.5,0.025,0.975))
    fitS[fitS$t==i & fitS$sero==s,3] <- data$Ncases[i,s]
  }
}
for(i in 1:data$Nt){
  fitTot[fitTot$t==i,3:5] <- quantile(chains$predCasesTot[,i], c(0.5,0.025,0.975))
  fitTot[fitTot$t==i,2] <- data$NcasesTot[i]
}

# Plot time-series fit
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

# Bt PLot
BtS <- data.frame(t=seq(1:data$Nt),med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt){
  for(s in 1:4){
    BtS[BtS$t==i,2:4] <- quantile(chains$B[,i], c(0.5,0.025,0.975))
  }
}
BtPlot <- ggplot(BtS, aes(t,med))+ geom_point(col='blue')+
  theme_minimal()+ ylim(0,NA) + theme(text=element_text(size=16))+
  theme(legend.position='none')+ ylab('Bt')+ xlab('Time')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU),fill='dodgerblue',alpha=0.5)
BtPlot


# Save plots
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
plot(plotfit1_5)
dev.off()
png(filename='6.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_1)
dev.off()
png(filename='7.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_2)
dev.off()
png(filename='8.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_3)
dev.off()
png(filename='9.png', width=20, height=20, res=400, units='cm')
plot(plotfit2_4)
dev.off()



