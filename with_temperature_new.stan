//----- 4-Serotype SEIR model  -----//
  // assumes immunity after 2 infections
// no age structure

data {
  
  int<lower=0> Nt; // N data time points
  int Ncases[Nt,4]; // Serotyped cases
  int NcasesTot[Nt]; // Non-serotyped cases
  int pop; // population
  real gamma; // birth/death rate
  vector[Nt] temp; //temperature
  real sigma; // 1/infectious period
  real omega; // 1/latent period
  real kappa; // human-mosquito-human transmission time
  vector[12] alpha; // dirichlet prior
  vector[4] pSERO; // serotype porportions  
  real phi; // ADE factor
  real zeta; // 1 / duration of the TCI period
  real delta; // relative reporting rate for primary infections
  real rho; // reporting rate
  real m; // briere function parameters
  real a; 
  real Tmin;
  real Tmax;
  
}


parameters {
  
  simplex[12] init; // initial conditions
  
}


transformed parameters {
  
  vector[Nt] casetot; // total reported cases
  matrix[Nt,4] cases; // serotype cases
  matrix[Nt-1,4] lam; // FOI
  matrix[Nt,4] Rt; // effective reproductive number
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
    real initV1 = init[3];
    vector[4] initMoP = init[4:7];
    vector[4] initMoS = init[8:11];
    real initV2 = init[12];
    real initV = initV1 + initV2 + sum(initMoP);
    
    
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
    e1[1, ] = to_row_vector((sigma/(sigma+omega))*initV1*pSERO); 
    e2[1, ] = to_row_vector((sigma/(sigma+omega))*initV2*pSERO);
    i1[1, ] = to_row_vector((omega/(sigma+omega))*initV1*pSERO);
    i2[1, ] = to_row_vector((omega/(sigma+omega))*initV2*pSERO);
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
    

    for(t in 1: Nt) cases[t,] = rho*omega*(e2[t,] + delta*e1[t,])*pop;
    for(t in 1: Nt) casetot[t] = sum(cases[t,]);
    for(t in 1: Nt) Rt[t,] = (su[t] + (sum(moS[t,])) - moS[t,])*B[t]/(sigma+kappa);
    
  
    
  }
}


model {

  // likelihood
  target += poisson_lpmf(Ncases[,1] | cases[,1]) + poisson_lpmf(Ncases[,2] | cases[,2]) + poisson_lpmf(Ncases[,3] | cases[,3]) + poisson_lpmf(Ncases[,4] | cases[,4]); 
  
  // Priors
  init ~ dirichlet(alpha);
    
  }
  
  
  generated quantities {
    
    real predCases[Nt,4];
    real predCasesTot[Nt];
    
    for(s in 1:4) predCases[,s] = poisson_rng(cases[,s]);
    predCasesTot = poisson_rng(casetot);
    
  }