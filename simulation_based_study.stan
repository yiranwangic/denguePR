//----- 4-Serotype SEIR model  -----//
  // assumes immunity after 2 infections
// no age structure

data {
  
  int<lower=0> Nt; // N data time points
  int Ncases[Nt,4]; // Serotyped cases
  int NcasesTot[Nt]; // Non-serotyped cases
  int pop; // population
  real gamma; // birth/death rate
  real temp[Nt]; //temperature
  real sigma; // 1/infectious period
  real omega; // 1/latent period
  real kappa; // human-mosquito-human transmission time
  vector[14] alpha; // dirichlet prior

  
}


parameters {
  
  real <upper=0> log_rho; // reporting rate
  simplex[14] init; // initial conditions
  real <upper=0> log_delta; // relative reporting rate for primary infection
  real <lower=0> a;
  real <lower=0> Tmin;//lower than this value, bt=0
  real <lower=0> Tmax; //higher than this value,bt=0
  real <upper=0> log_zeta; //TCI Period
  real <lower=1> phi; //ADE factor
  
  
}

transformed parameters {
  
  real rho = exp(log_rho);
  real delta = exp(log_delta);
  real zeta =exp(log_zeta);
  real initS = init[1];
  real initR = init[2];
  vector[4] initV1 = init[3:6];
  vector[4] initMoS = init[7:10];
  vector[4] initV2 = init[11:14];
  real initV = sum(initV1) + sum(initV2);
  vector[Nt] casetot; // total reported cases
  matrix[Nt,4] cases; // serotype cases
  matrix[Nt,4] lam; // FOI
  matrix[Nt,4] mz; // Mosquito
  real lamtot[Nt];
  real Rt[Nt,4];
  vector[Nt] su;
  vector[Nt] B;
  for(t in 1:Nt) B[t]= a*(temp[t]^2)*((temp[t]-Tmin)^2)*(Tmax-temp[t]);
  
  
  
  
  
  {  
    
    
    
    matrix[Nt,4] e1;
    matrix[Nt,4] i1;
    matrix[Nt,4] e2;
    matrix[Nt,4] i2;
    matrix[Nt,4] moP;
    matrix[Nt,4] moS;
    vector[Nt] r;
    
    
    
    
    
    
    // initial conditions
    
    su[1] = initS;
    r[1] = initR;
    
    
    for(s in 1:4){
      moP[1,s] = (omega*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[s];
      moS[1,s] = initMoS[s];
       e1[1,s] = (zeta*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[s];
       e2[1,s] = (sigma/(sigma+omega))*initV2[s];
       i1[1,s] = (zeta*omega/(omega*sigma+zeta*sigma+zeta*omega))*initV1[s];
       i2[1,s] = (omega/(sigma+omega))*initV2[s];
       mz[1,s] = B[1]*(i1[1,s]+phi*i2[1,s]);
       
    }
    
    
    
    // loop through time & age
    for(t in 1:(Nt-1)){
      
      for(s in 1:4) mz[t+1,s] = mz[t,s]+ B[t]*(i1[t,s] + phi*i2[t,s]) - kappa*mz[t,s];
      for(s in 1:4) lam[t,s] = kappa*mz[t,s];
      lamtot[t] = sum(lam[t,]);
      
      su[t+1] = su[t] - su[t]*sum(lam[t,]) - gamma*su[t] + gamma;
      for(s in 1:4) e1[t+1,s] = e1[t,s] + su[t]*lam[t,s] - omega*e1[t,s] - gamma*e1[t,s];
      for(s in 1:4) i1[t+1,s] = i1[t,s] + omega*e1[t,s] - sigma*i1[t,s] - gamma*i1[t,s];
      for(s in 1:4) moP[t+1,s] = moP[t,s] + sigma*i1[t,s] - moP[t,s]*zeta - gamma*moP[t,s];
      for(s in 1:4) moS[t+1,s] = moS[t,s] + moP[t,s]*zeta - moS[t,s]*(lamtot[t]-lam[t,s]) - gamma*moS[t,s];
      for(s in 1:4) e2[t+1,s] = e2[t,s] + (sum(moS[t,])-moS[t,s])*lam[t,s] - omega*e2[t,s] - gamma*e2[t,s];
      for(s in 1:4) i2[t+1,s] = i2[t,s] + omega*e2[t,s] - sigma*i2[t,s] - gamma*i2[t,s];
      r[t+1] = r[t] + sigma*sum(i2[t,]) - gamma*r[t];
    }  
    for(t in 1: Nt) for(s in 1:4) cases[t,s] = rho*omega*(e2[t,s] + delta*e1[t,s])*pop;
    for(t in 1: Nt) casetot[t] = sum(cases[t,]);
    
  }
}


model {
  
  // Priors
  init ~ dirichlet(alpha);
  log_rho ~ normal(-1.6,0.5);
  log_delta ~ normal(-3,0.25);
  a ~ normal(4e-7,5e-8);
  Tmin ~ normal(13,2);
  Tmax ~ normal(37,2);
  log_zeta ~ normal(-6,1);
  phi ~ normal(1.5,0.5);
  
  // likelihood
  NcasesTot ~ poisson(casetot); // neg-binomial total cases
  for(s in 1:4) Ncases[,s] ~ poisson(cases[,s]);
  
}


generated quantities {
  
  real predCases[Nt,4];
  real predCasesTot[Nt];
  
  for(s in 1:4) predCases[,s] = poisson_rng(cases[,s]);
  predCasesTot = poisson_rng(casetot);
  
}