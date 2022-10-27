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
  real sigma; // 1/infectious period
  real omega; // 1/latent period
  real kappa; // human-mosquito-human transmission time
  real phi; //ADE factor
  real br[Nt2]; // birth rate
  real er[Nt2]; // exit rate
  int K; // N compartments
  vector[K] alpha; // dirichlet prior
  int indL[Nt]; // indices for time step
  int indU[Nt];
  int ind[Nt2];

}


parameters {
  
  real <upper=0> log_rho; // reporting rate
  simplex[K] init; // initial conditions
  real log_B[Nt]; // transmission rate
  real <lower=0> log_phiNB; //variance for beta
  real <lower=0> delta; //relative reporting rate for primary infections
}


transformed parameters {
  
  real rho = exp(log_rho);
  real initS = init[1];
  real initR = init[2];
  vector[4] initMo = init[3:6];
  vector[4] initE1 = init[7:10];
  vector[4] initI1 = init[11:14];
  vector[4] initE2 = init[15:18];
  vector[4] initI2 = init[19:22];
  
  vector[Nt] Ctot; // total reported cases
  real B[Nt] = exp(log_B); // transmission rate
  matrix[Nt,4] C; // total cases by serotype
  matrix[Nt,4] rC; // reported serotype cases
  matrix[Nt2-1,4] lam; // FOI
  matrix[Nt2,4] Mz; // Mosquito
  simplex[4] pSero[Nt];
  real Rt[Nt,4];
  vector[Nt2] S;
  
  {
    
    matrix[Nt2,4] E1;
    matrix[Nt2,4] I1;
    matrix[Nt2,4] E2;
    matrix[Nt2,4] I2;
    matrix[Nt2,4] Mo;
    vector[Nt2] R;
    real lamtot[Nt2];
    
    
    // initial conditions
    S[1] = init[1];
    R[1] = init[2];
    for(s in 1:4){
      Mo[1,s] = initMo[s];
      E1[1,s] = initE1[s];
      E2[1,s] = initE2[s];
      I1[1,s] = initI1[s];
      I2[1,s] = initI2[s];
      Mz[1,s] = B[1]*(initI1[s]+phi*initI2[s]);
    }
    
    
    // loop through time & age
    for(t in 1:(Nt2-1)){
      
      for(s in 1:4) Mz[t+1,s] = Mz[t,s] + B[ind[t]]*(I1[t,s] + phi*I2[t,s]) - kappa*Mz[t,s];
      for(s in 1:4) lam[t,s] = kappa*Mz[t,s];
      lamtot[t] = sum(lam[t,]);
      
      S[t+1] = S[t] - S[t]*sum(lam[t,]) - er[t]*S[t] + br[t];
      for(s in 1:4) E1[t+1,s] = E1[t,s] + S[t]*lam[t,s] - omega*E1[t,s] - er[t]*E1[t,s];
      for(s in 1:4) I1[t+1,s] = I1[t,s] + omega*E1[t,s] - sigma*I1[t,s] - er[t]*I1[t,s];
      for(s in 1:4) Mo[t+1,s] = Mo[t,s] + sigma*I1[t,s] - Mo[t,s]*(lamtot[t]-lam[t,s]) - er[t]*Mo[t,s];
      for(s in 1:4) E2[t+1,s] = E2[t,s] + (sum(Mo[t,]) - Mo[t,s])*lam[t,s] - omega*E2[t,s] - er[t]*E2[t,s];
      for(s in 1:4) I2[t+1,s] = I2[t,s] + omega*E2[t,s] - sigma*I2[t,s] - er[t]*I2[t,s];
      R[t+1] = R[t] + sigma*sum(I2[t,]) - er[t]*R[t];
    }
    
    for(s in 1:4) for(t in 1:Nt) C[t,s] = rho*omega*(sum(E2[indL[t]:indU[t],s]) + delta*sum(E1[indL[t]:indU[t],s]))*pop[t] +0.001;
    for(t in 1:Nt) Ctot[t] = sum(C[t,]); 
    for(t in 1:Nt) for(s in 1:4) rC[t,s] = C[t,s]*pSample[t] + 0.001;
    for(t in 1:Nt) for(s in 1:4) pSero[t,s] = sum(I1[indL[t]:indU[t],s]+I2[indL[t]:indU[t],s])/((sum(I1[indL[t]:indU[t],1:4])+sum(I2[indL[t]:indU[t],1:4])));
    for(t in 1:Nt) for(s in 1:4) Rt[t,s] = (S[t]+sum(Mo[t,])-Mo[t,s])*B[t]/sigma;
    
  }
  
}



model {
  
  // Priors
  log_B[1] ~ normal(-0.5,0.7);
  init ~ dirichlet(alpha);
  log_rho ~ normal(-2,0.5);
  delta ~ normal(0.3,0.1);
  log_phiNB ~ normal(0.5,0.5);
    
  // likelihood
  for(t in 2:Nt) log_B[t] ~ normal(log_B[t-1], 0.05);
  NcasesTot ~ neg_binomial_2(Ctot, exp(log_phiNB)); // neg-binomial total cases
  for(s in 1:4) Ncases[,s] ~ neg_binomial_2(rC[,s], exp(log_phiNB));

}


generated quantities {
  
  real predCases[Nt,4];
  real predCasesTot[Nt];
  
  for(s in 1:4) predCases[,s] = neg_binomial_2_rng(rC[,s], exp(log_phiNB));
  predCasesTot = neg_binomial_2_rng(Ctot, exp(log_phiNB));

}