#---- Fit SEIR model to San Juan serotype case data -----#

# libraries
rm(list=ls())
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

# inputs for model
data <- list()
data$Ncases <- Ncases
data$NcasesTot <- df$total_cases
data$Nt <- nrow(data$Ncases)
data$Nt2 <- data$Nt*7
# population for each year of the time-series
data$pop <- c(rep(2327000,35),
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
data$phi <- 1.5
data$zeta <- 1/0.5/365
# annual average birth rate for each year of the time-series
data$br<- c(rep(4.96e-05,35*7),
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
data$er<- c(rep(2.84e-05,35*7),
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
data$alpha <- c(20,60,rep(5.2,3),rep(5,3),rep(0.2,3))/1
indL <- seq(1,data$Nt*7,7)
indU <- seq(7,data$Nt*7+6,7)
data$indL <- indL
data$indU <- indU
data$ind <- sort(rep(seq(1,data$Nt),7))
data$pSample <- rowSums(data$Ncases)/data$NcasesTot
data$pSample[is.na(data$pSample)] <- 0
data$temp <- df$TAVG


# chain starting values
alpha <- c(20,60,rep(5.2,3),rep(5,3),rep(0.2,3))/1
initI <- rdirichlet(100,alpha)
ii <- as.vector(initI[1,])
inits_chain1 <- list(log_rho=-1.6, init=ii, log_delta=-3,a=1.5e-7,Tmin=12, Tmax=36, log_est=-14)  
inits_chain2 <- list(log_rho=-1.7, init=ii, log_delta=-2, a=2e-7, Tmin=13, Tmax=37, log_est=-13.5)
inits_chain3 <- list(log_rho=-1.8, init=ii, log_delta=-2.5,a=3e-7,Tmin=14, Tmax=38, log_est=-14.5)

# fit model
check_cmdstan_toolchain(fix=T)
set_cmdstan_path()
mod <- cmdstan_model('fit_to_san_juan_data.stan', pedantic=T)
fit <- mod$sample(data=data, chains=3, parallel_chains=3, iter_sampling=900, refresh=10, iter_warmup=100, save_warmup=TRUE, init=list(inits_chain1,inits_chain2,inits_chain3))
stanfit <- rstan::read_stan_csv(fit$output_files())

# Check convergence
chains <- rstan::extract(stanfit)
traceplot1 <- traceplot(stanfit, pars=c('init','lp__'),inc_warmup=T)
traceplot2 <- traceplot(stanfit, pars=c('rho','delta','a','Tmin','Tmax','log_est'),inc_warmup=T)

summary(stanfit,pars=c('init','lp__'),probs = c(0.025, 0.975))$summary
summary(stanfit,pars=c('rho','delta','a','Tmin','Tmax','log_est'),probs = c(0.025, 0.975))$summary

# plot time-series fit
fitS <- data.frame(date=as.Date(dates),t=rep(seq(1:data$Nt),4),sero=sort(rep(seq(1,4),data$Nt)),obs=NA, med=NA,ciL=NA,ciU=NA)
fitTot <- data.frame(date=as.Date(dates),t=seq(1,data$Nt),obs=data$NcasesTot, med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt){
  for(s in 1:4){
    fitS[fitS$t==i & fitS$sero==s,5:7] <- quantile(chains$predCases[,i,s], c(0.5,0.025,0.975))
    fitS[fitS$t==i & fitS$sero==s,4] <- data$Ncases[i,s]
  }
  fitTot[i,4:6] <- quantile(chains$predCasesTot[,i], c(0.5,0.025,0.975))
}

# Plots
fitPlotS <- ggplot(fitS, aes(date,obs))+ geom_point(alpha=0.4, col='grey40')+
  theme_minimal()+ theme(legend.position='none')+ theme(text=element_text(size=16))+
  geom_line(aes(date,med,col=factor(sero)))+ geom_ribbon(aes(date, ymin=ciL,ymax=ciU,fill=factor(sero)),alpha=0.3)+
  xlab('Time')+ ylab('Cases')+ theme(legend.position='none')+ facet_wrap(~sero,scale='free_y')
fitPlotS
fitPlotTot <- ggplot(fitTot, aes(date,obs))+ geom_point(alpha=0.4, col='grey40')+
  theme_minimal()+ geom_line(aes(date,med), col='blue')+ xlab('Time')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU), fill='dodgerblue', alpha=0.5)+ ylab('Cases')+
  theme(text=element_text(size=16))
fitPlotTot

# Epidemiological parameter estimates
BtS <- data.frame(date=as.Date(dates),t=seq(1:data$Nt),med=NA,ciL=NA,ciU=NA)
Rt <- data.frame(date=as.Date(dates),t=rep(seq(1:data$Nt),4),sero=sort(rep(seq(1,4),data$Nt)),med=NA,ciL=NA,ciU=NA)
lam <- data.frame(date=seq(1:data$Nt2),t=rep(seq(1:data$Nt2),4),sero=sort(rep(seq(1,4),data$Nt)),med=NA,ciL=NA,ciU=NA)
for(i in 1:data$Nt){
  for(s in 1:4){
    BtS[BtS$t==i,3:5] <- quantile(chains$B[,i], c(0.5,0.025,0.975))
    Rt[Rt$t==i & Rt$sero==s,4:6] <- quantile(chains$Rt[,i,s], c(0.5,0.025,0.975))
    
  }
}
for(i in 1:data$Nt2) for(s in 1:4){
  lam[lam$t==i & lam$sero==s,4:6] <- quantile(chains$lam[,i,s], c(0.5,0.025,0.975))
}
BtPlot <- ggplot(BtS, aes(date,med))+ geom_point(col='blue')+
  theme_minimal()+ ylim(0,NA) + theme(text=element_text(size=16))+
  theme(legend.position='none')+ ylab('Bt')+ xlab('Time')+
  geom_ribbon(aes(ymin=ciL,ymax=ciU),fill='dodgerblue',alpha=0.5)
BtPlot
RtPlot <- ggplot(Rt, aes(date,med,col=factor(sero)))+ geom_point()+
  facet_wrap(~sero)+ theme_minimal()+ ylim(0,NA)+ theme(text=element_text(size=16))+
  geom_linerange(aes(ymin=ciL,ymax=ciU,fill=factor(sero)), alpha=0.3)+
  theme(legend.position='none')+ ylab('Rt')+ xlab('Time')+
  geom_hline(yintercept=1, linetype='dashed')
RtPlot
lamPlot <- ggplot(lam, aes(date,med))+ geom_line(aes(col=factor(sero)))+
  facet_wrap(~sero, scales='free_y')+ theme_minimal()+ ylim(0,NA)+ theme(text=element_text(size=16))+
  geom_ribbon(aes(ymin=ciL,ymax=ciU,fill=factor(sero)), alpha=0.3)+
  theme(legend.position='none')+ ylab('FOI')+ xlab('days')
lamPlot

# Save results
setwd('/Users/wangyr/Desktop/results')
saveRDS(data,'Data.RDS')
saveRDS(chains, 'Chains.RDS')
png(filename='FitTotalCases.png', width=20, height=14, res=400, units='cm')
plot(fitPlotTot)
dev.off()
png(filename='FitSeroCases.png', width=20, height=14, res=400, units='cm')
plot(fitPlotS)
dev.off()
png(filename='RtPlot.png', width=20, height=15, res=400, units='cm')
plot(RtPlot)
dev.off()
png(filename='BtPlot.png', width=20, height=15, res=400, units='cm')
plot(BtPlot)
dev.off()
png(filename='FOIPlot.png', width=20, height=15, res=400, units='cm')
plot(lamPlot)
dev.off()
png(filename='traceplot1.png', width=30, height=14, res=400, units='cm')
plot(traceplot1)
dev.off()
png(filename='traceplot2.png', width=20, height=14, res=400, units='cm')
plot(traceplot2)
dev.off()