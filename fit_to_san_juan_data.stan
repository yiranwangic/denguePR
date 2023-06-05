//----- 4-Serotype SEIR model  -----//
  // assumes immunity after 2 infections
// no age structure

data {
  
  int<lower=0> Nt; // N data time points
  int Nt2; // N time points for modelling
  int Ncases[Nt,4]; // Serotyped cases
  int NcasesTot[Nt]; // Non-serotyped cases
  real pSample[Nt]; // proportion of cases serotyped
  int pop[Nt]; // population
  real temp[Nt]; //temperature
  real sigma; // 1/infectious period
  real omega; // 1/latent period
  real kappa; // human-mosquito-human transmission time
  real br[Nt2]; // birth rate
  real er[Nt2]; // exit rate
  vector[11] alpha; // dirichlet prior
  int indL[Nt]; // indices for time step
  int indU[Nt];
  int ind[Nt2];
  real phi; // ADE factor
  real zeta; // TCI

  
}


parameters {
  
  real <upper=0> log_rho; // reporting rate
  simplex[11] init; // initial conditions
  real <upper=0> log_delta; // relative reporting rate for primary infection
  real <lower=0> a;
  real <lower=0> Tmin;//lower than this value, bt=0
  real <lower=0> Tmax; //higher than this value,bt=0
  real <upper=0> log_est;

  
}

transformed parameters {
  
  real rho = exp(log_rho);
  real delta = exp(log_delta);
  real est = exp(log_est); // reseed DENV-3
  real initS = init[1];
  real initR = init[2];
  vector[3] initV1 = init[3:5];
  vector[3] initMoS = init[6:8];
  vector[3] initV2 = init[9:11];
  real initV = sum(initV1) + sum(initV2);
  vector[Nt] Ctot; // total reported cases
  matrix[Nt,4] C; // serotype cases
  matrix[Nt,4] rC; // reported serotype cases
  matrix[Nt2-1,4] lam; // FOI
  matrix[Nt2,4] Mz; // Mosquito
  simplex[4] pSero[Nt];
  real Rt[Nt,4];
  real lamtot[Nt2];
  vector[Nt2] S;
  vector[Nt] Bt; // temperature-dependent bt
  for(t in 1:1196) Bt[t]= a*(temp[t]^2)*((temp[t]-Tmin)^2)*(Tmax-temp[t]);
  
  
  
  
  
  {  
    
    
    
    matrix[Nt2,4] E1;
    matrix[Nt2,4] I1;
    matrix[Nt2,4] E2;
    matrix[Nt2,4] I2;
    matrix[Nt2,4] MoP;
    matrix[Nt2,4] MoS;
    vector[Nt2] R;
    
    
    
    
    
    
    // initial conditions
    
    S[1] = init[1];
    R[1] = init[2];
    MoP[1,1] = (omega*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[1];
    MoP[1,2] = (omega*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[2];
    MoP[1,3] = 0;
    MoP[1,4] = (omega*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[3];
    MoS[1,1] = initMoS[1];
    MoS[1,2] = initMoS[2];
    MoS[1,3] = 0;
    MoS[1,4] = initMoS[3];
    E1[1,1] = (zeta*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[1];
    E1[1,2] = (zeta*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[2];
    E1[1,3] = 0;
    E1[1,4] = (zeta*sigma/(omega*sigma+zeta*sigma+zeta*omega))*initV1[3];
    E2[1,1] = (sigma/(sigma+omega))*initV2[1];
    E2[1,2] = (sigma/(sigma+omega))*initV2[2];
    E2[1,3] = 0;
    E2[1,4] = (sigma/(sigma+omega))*initV2[3];
    I1[1,1] = (zeta*omega/(omega*sigma+zeta*sigma+zeta*omega))*initV1[1];
    I1[1,2] = (zeta*omega/(omega*sigma+zeta*sigma+zeta*omega))*initV1[2];
    I1[1,3] = 0;
    I1[1,4] = (zeta*omega/(omega*sigma+zeta*sigma+zeta*omega))*initV1[3];
    I2[1,1] = (omega/(sigma+omega))*initV2[1];
    I2[1,2] = (omega/(sigma+omega))*initV2[2];
    I2[1,3] = 0;
    I2[1,4] = (omega/(sigma+omega))*initV2[3];
    
    for(s in 1:4){
      Mz[1,s] = Bt[1]*(I1[1,s]+phi*I2[1,s]);
    }
    
    
    
    // loop through time & age
    for(t in 1:(Nt2-1)){
      
      for(s in 1:4) Mz[t+1,s] = Mz[t,s] + Bt[ind[t]]*(I1[t,s] + phi*I2[t,s]) - kappa*Mz[t,s];
      for(s in 1:4) lam[t,s] = kappa*Mz[t,s];
      lamtot[t] = sum(lam[t,]);
      
      S[t+1] = S[t] - S[t]*sum(lam[t,]) - er[t]*S[t] + br[t];
      for(s in 1:4) E1[t+1,s] = E1[t,s] + S[t]*lam[t,s] - omega*E1[t,s] - er[t]*E1[t,s];
      if(t==398*7) I1[t,3] = I1[t,3] + est; //reseed DENV-3
      for(s in 1:4) I1[t+1,s] = I1[t,s] + omega*E1[t,s] - sigma*I1[t,s] - er[t]*I1[t,s];
      for(s in 1:4) MoP[t+1,s] = MoP[t,s] + sigma*I1[t,s] - zeta*MoP[t,s] - er[t]*MoP[t,s];
      for(s in 1:4) MoS[t+1,s] = MoS[t,s] + zeta*MoP[t,s] - MoS[t,s]*(lamtot[t]-lam[t,s]) - er[t]*MoS[t,s];
      for(s in 1:4) E2[t+1,s] = E2[t,s] + (sum(MoS[t,]) - MoS[t,s])*lam[t,s] - omega*E2[t,s] - er[t]*E2[t,s];
      for(s in 1:4) I2[t+1,s] = I2[t,s] + omega*E2[t,s] - sigma*I2[t,s] - er[t]*I2[t,s];
      R[t+1] = R[t] + sigma*sum(I2[t,]) - er[t]*R[t];
    }
    
    for(s in 1:4) for(t in 1:Nt) C[t,s] = rho*omega*(sum(E2[indL[t]:indU[t],s]) + delta*sum(E1[indL[t]:indU[t],s]))*pop[t] +0.001;
    for(t in 1:Nt) Ctot[t] = sum(C[t,]); 
    for(t in 1:Nt) for(s in 1:4) rC[t,s] = C[t,s]*pSample[t] + 0.001;
    for(t in 1:Nt) for(s in 1:4) pSero[t,s] = sum(I1[indL[t]:indU[t],s]+I2[indL[t]:indU[t],s])/((sum(I1[indL[t]:indU[t],1:4])+sum(I2[indL[t]:indU[t],1:4])));
    for(t in 1:Nt) for(s in 1:4) Rt[t,s] = (S[t]+sum(MoS[t,])-MoS[t,s])*Bt[t]/sigma;
  }
  
}



model {
  
  // Priors
  init ~ dirichlet(alpha);
  log_rho ~ normal(-2.5,0.4);
  log_delta ~ normal(-2,0.25);
  a ~ normal(1.5e-7,5e-8);
  Tmin ~ normal(13,2);
  Tmax ~ normal(37,2);
  log_est ~ normal(-14,0.5); // reseed DENV-3

  
  // likelihood
  target+= poisson_lpmf(Ncases[,1] | rC[,1]) + poisson_lpmf(Ncases[,2] | rC[,2]) + poisson_lpmf(Ncases[,3] | rC[,3]) + poisson_lpmf(Ncases[,4] | rC[,4]);

}


generated quantities {
  
  real predCases[Nt,4];
  real predCasesTot[Nt];


  for(s in 1:4) predCases[,s] = poisson_rng(rC[,s]);
  predCasesTot = poisson_rng(Ctot);
  
}