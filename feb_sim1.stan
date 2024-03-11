//----- 4-Serotype SEIR model  -----//
  // assumes immunity after 2 infections
// no age structure

data {
  
  int<lower=0> Nt; // N data time points
  int<lower=0> Nt2; // N data time points
  int indL[Nt2]; // indices for daily -> weekly
  int indU[Nt2];
  int Ncases[Nt2,4]; // Serotyped cases
  int NcasesTot[Nt2]; // Non-serotyped cases
  int pop; // population
  real gamma; // birth/death rate
  vector[Nt] temp; //temperature
  real sigma; // 1/infectious period
  real omega; // 1/latent period
  real kappa; // human-mosquito-human transmission time
  vector[14] alpha; // dirichlet prior
  real phi; // ADE factor
  real zeta; // 1 / duration of the TCI period
  real delta; // relative reporting rate for primary infections
  real rho; // reporting rate
  real m;
  
  
}


parameters {
  
  simplex[14] init; // initial conditions
  real <upper=0> log_a;
  real <lower=0> Tmin;
  real <lower=0> Tmax;
  real <lower=0, upper=1> pv1_s1;
  real <lower=0, upper=1> pv1_s2;
  real <lower=0, upper=1> pv1_s3;
  real <lower=0, upper=1> pv1_s4;
  
}


transformed parameters {
  
  real a = exp(log_a);
  vector[Nt2] casetot; // total reported cases
  matrix[Nt2,4] cases; // serotype cases
  matrix[Nt-1,4] lam; // FOI
  vector[Nt] B; // temperature-dependent beta_t
  
  // General briere function
  for (t in 1:Nt) {
    if (temp[t] > Tmax || temp[t] < Tmin) {
      B[t] = 0;
    } else {
      B[t] = a * temp[t] * (temp[t] - Tmin) * ((Tmax - temp[t])^(1 / m));
    }
  }
  
  
  { // this { defines a block within which variables declared are local (can't have lower or upper)
    real initS = init[1];
    real initR = init[2];
    real initV_s1 = init[3];
    real initV_s2 = init[4];
    real initV_s3 = init[5];
    real initV_s4 = init[6];
    vector[4] initMoP = init[7:10];
    vector[4] initMoS = init[11:14];
    
    
    vector[Nt] su;  
    matrix[Nt,4] e1;
    matrix[Nt,4] i1;
    matrix[Nt,4] e2;
    matrix[Nt,4] i2;
    matrix[Nt,4] moP;
    matrix[Nt,4] moS;
    vector[Nt] r;
    matrix[Nt,4] mz; // Mosquito compartment accounting for the EIP
    real lamtot[Nt-1];
    

// initial conditions
    
    su[1] = initS;
    r[1] = initR;
    

    moP[1, ] = to_row_vector(initMoP);
    moS[1, ] = to_row_vector(initMoS);
    
    e1[1, 1] = (sigma/(sigma+omega))*initV_s1*pv1_s1; 
    i1[1, 1] = (omega/(sigma+omega))*initV_s1*pv1_s1;
    
    e1[1, 2] = (sigma/(sigma+omega))*initV_s2*pv1_s2; 
    i1[1, 2] = (omega/(sigma+omega))*initV_s2*pv1_s2;
    
    e1[1, 3] = (sigma/(sigma+omega))*initV_s3*pv1_s3; 
    i1[1, 3] = (omega/(sigma+omega))*initV_s3*pv1_s3;
    
    e1[1, 4] = (sigma/(sigma+omega))*initV_s4*pv1_s4; 
    i1[1, 4] = (omega/(sigma+omega))*initV_s4*pv1_s4;
    
    e2[1, 1] = (sigma/(sigma+omega))*initV_s1*(1-pv1_s1);
    i2[1, 1] = (omega/(sigma+omega))*initV_s1*(1-pv1_s1);
    
    e2[1, 2] = (sigma/(sigma+omega))*initV_s2*(1-pv1_s2);
    i2[1, 2] = (omega/(sigma+omega))*initV_s2*(1-pv1_s2);
    
    e2[1, 3] = (sigma/(sigma+omega))*initV_s3*(1-pv1_s3);
    i2[1, 3] = (omega/(sigma+omega))*initV_s3*(1-pv1_s3);
    
    e2[1, 4] = (sigma/(sigma+omega))*initV_s4*(1-pv1_s4);
    i2[1, 4] = (omega/(sigma+omega))*initV_s4*(1-pv1_s4);
    
    mz[1, ] = B[1] * ((i1[1, ] + phi ) .* (i2[1, ]));
    
    
    
    // loop through time & age
    for(t in 1:(Nt-1)){
      
      mz[t+1, ] = mz[t, ] + (B[t] * (i1[t, ] + phi*i2[t, ])) - (kappa*mz[t, ]);

      lam[t, ] = kappa*mz[t, ]; 
      
      lamtot[t] = sum(lam[t,]);
      
      su[t+1] = su[t] - su[t]*sum(lam[t,]) - gamma*su[t] + gamma;

      e1[t+1, ] = e1[t, ] + su[t]*lam[t, ] - omega*e1[t, ] - gamma*e1[t, ];

      i1[t+1, ] = i1[t, ] + omega*e1[t, ] - sigma*i1[t, ] - gamma*i1[t, ]; 

      moP[t+1, ] = moP[t, ] + sigma*i1[t, ] - moP[t, ]*zeta - gamma*moP[t, ]; 

      moS[t+1, ] = moS[t, ] + moP[t, ]*zeta - moS[t, ].*(lamtot[t] - lam[t, ]) - gamma*moS[t, ]; 

      e2[t+1, ] = e2[t, ] + ((sum(moS[t, ]) - moS[t, ]) .* lam[t, ]) - omega*e2[t, ] - gamma*e2[t, ]; 

      i2[t+1, ] = i2[t, ] + omega*e2[t, ] - sigma*i2[t, ] - gamma*i2[t, ]; 

      r[t+1] = r[t] + sigma*sum(i2[t,]) - gamma*r[t];
    }  
    

    for(t in 1: Nt2) for(s in 1:4) cases[t,s] = rho*omega*(sum(e2[indL[t]:indU[t],s]) + delta*sum(e1[indL[t]:indU[t],s]))*pop;
    for(t in 1: Nt2) casetot[t] = sum(cases[t,]);

  
    
  }
}


model {

  // likelihood
  target += poisson_lpmf(Ncases[,1] | cases[,1]) + poisson_lpmf(Ncases[,2] | cases[,2]) + poisson_lpmf(Ncases[,3] | cases[,3]) + poisson_lpmf(Ncases[,4] | cases[,4]); 
  
  // Priors
  init ~ dirichlet(alpha);
  log_a ~ normal(-6.5,0.5);
  Tmin ~ normal(17.8,3.5);
  Tmax ~ normal(34.6,1);
  pv1_s1 ~ normal(0.5, 0.2);
  pv1_s2 ~ normal(0.5, 0.2);
  pv1_s3 ~ normal(0.5, 0.2);
  pv1_s4 ~ normal(0.5, 0.2);
  
  }
  
  
  generated quantities {
    
    real predCases[Nt2,4];
    real predCasesTot[Nt2];
    
    for(s in 1:4) predCases[,s] = poisson_rng(cases[,s]);
    predCasesTot = poisson_rng(casetot);
    
  }