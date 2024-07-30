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
  real temp[Nt2]; //temperature
  real sigma; // 1/infectious period
  real omega; // 1/latent period
  real kappa; // human-mosquito-human transmission time
  real br[Nt2]; // birth rate
  real er[Nt2]; // exit rate
  vector[14] alpha1; // dirichlet prior
  int indL[Nt]; // indices for time step
  int indU[Nt];
  real m;
  
}


parameters {
  
  real <upper=0> log_rho; // reporting rate
  real <upper=0> log_delta; // relative reporting rate for primary infection
  simplex[14] init; // initial conditions
  real <lower=0> a;
  real <lower=0> Tmin;//lower than this value, bt=0
  real <lower=0> Tmax; //higher than this value,bt=0
  real <upper=0> log_est;
  real <lower=0, upper=1> pv1;
  real <lower=1> phi;
  real <upper=0> log_zeta;
  
}

transformed parameters {
  
  real zeta = exp(log_zeta);
  real rho = exp(log_rho);
  real delta = exp(log_delta);
  real est = exp(log_est); // reseed DENV-3

  vector[Nt] Ctot; // total reported cases
  matrix[Nt,4] C; // serotype cases
  matrix[Nt,4] rC; // reported serotype cases
  matrix[Nt2-1,4] lam; // FOI
  real lamtot[Nt2-1];

  vector[Nt2] Bt; // temperature-dependent bt
  
  // General briere function
  for (t in 1:Nt2) {
    if (temp[t] > Tmax || temp[t] < Tmin) {
      Bt[t] = 0;
    } else {
      Bt[t] = a * temp[t] * (temp[t] - Tmin) * ((Tmax - temp[t])^(1 / m));
    }
  }


  
  {  
    
   real initS = init[1];
   real initR = init[2];
   vector[4] initV = init[3:6];
   vector[4] initMoS = init[7:10];
   vector[4] initMoP = init[11:14];    
    
    vector[Nt2] S;
    matrix[Nt2,4] E1;
    matrix[Nt2,4] I1;
    matrix[Nt2,4] E2;
    matrix[Nt2,4] I2;
    matrix[Nt2,4] MoP;
    matrix[Nt2,4] MoS;
    vector[Nt2] R;
    
    matrix[Nt2,4] Mz; // Mosquito
    
    

    // initial conditions
    
    S[1] = init[1];
    R[1] = init[2];
    MoP[1,] = to_row_vector(initMoP);
    MoS[1,] = to_row_vector(initMoS);
    E1[1,] = (sigma/(sigma+omega))*pv1*to_row_vector(initV);
    E2[1,] = (sigma/(sigma+omega))*(1-pv1)*to_row_vector(initV);
    I1[1,] = (omega/(sigma+omega))*pv1*to_row_vector(initV);
    I2[1,] = (omega/(sigma+omega))*(1-pv1)*to_row_vector(initV);
    
    for(s in 1:4){
      Mz[1,s] = Bt[1]*(I1[1,s]+phi*I2[1,s]);
    }
    
    
    
    // loop through time & age
    for(t in 1:(Nt2-1)){
      
      for(s in 1:4) Mz[t+1,s] = Mz[t,s] + Bt[t]*(I1[t,s] + phi*I2[t,s]) - kappa*Mz[t,s];
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
    
    for(s in 1:4) for(t in 1:Nt) C[t,s] = omega*(rho*sum(E2[indL[t]:indU[t],s]) + delta*sum(E1[indL[t]:indU[t],s]))*pop[t] +0.001;
    for(t in 1:Nt) Ctot[t] = sum(C[t,]); 
    for(t in 1:Nt) for(s in 1:4) rC[t,s] = C[t,s]*pSample[t] + 0.001;
  }
  
}



model {
  
  // Priors
  init ~ dirichlet(alpha1);
  log_rho ~ normal(-2.3,0.3);
  log_delta ~ normal(-4,0.5);
  a ~ normal(1.5e-7,5e-8);
  Tmin ~ normal(17.8,3.5);
  Tmax ~ normal(34.6,1);
  log_est ~ normal(-14,0.5); // reseed DENV-3
  pv1 ~ normal(0.5,0.1);
  phi ~ normal(1.5,0.3);
  log_zeta ~ normal(-6.3,1);


  
  // likelihood
  target+= poisson_lpmf(Ncases[,1] | rC[,1]) + poisson_lpmf(Ncases[,2] | rC[,2]) + poisson_lpmf(Ncases[,3] | rC[,3]) + poisson_lpmf(Ncases[,4] | rC[,4]);

}


generated quantities {
  
  real predCases[Nt,4];
  real predCasesTot[Nt];


  for(s in 1:4) predCases[,s] = poisson_rng(rC[,s]);
  predCasesTot = poisson_rng(Ctot);
  
}