// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only


//iCAR function
functions {
  real icar_normal_lpdf(vector bb, int nroutes, int[] node1, int[] node2) {
    return -0.5 * dot_self(bb[node1] - bb[node2])
      + normal_lpdf(sum(bb) | 0, 0.001 * nroutes); //soft sum to zero constraint on phi
 }
}


data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;
 
  int<lower=0> count[ncounts];              // count observations
  int<lower=1> year[ncounts]; // year index
  int<lower=1> route[ncounts]; // route index
  int<lower=0> firstyr[ncounts]; // first year index
  int<lower=1> observer[ncounts];              // observer indicators

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  int<lower=1, upper=nroutes> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=nroutes> node2[N_edges];  // and node1[i] < node2[i]


}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
  real lambda[ncounts];             // Poisson means
  
  vector[nroutes] beta_raw;
  real BETA; 

  vector[nroutes] alpha_raw;
  real ALPHA; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts

  
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;

// covariate effect on intercepts and slopes
   beta = (sdbeta*beta_raw) + BETA;
   alpha = (sdalpha*alpha_raw) + ALPHA;
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  }
  
model {

  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ student_t(4,0,1); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of slope variation

  beta_raw ~ icar_normal_lpdf(nroutes, node1, node2);
  alpha_raw ~ icar_normal_lpdf(nroutes, node1, node2);

  //sum to zero constraints
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  sum(beta_raw) ~ normal(0,0.001*nroutes);
  
}

 generated quantities {

     vector[ncounts] log_lik;
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  

 }

