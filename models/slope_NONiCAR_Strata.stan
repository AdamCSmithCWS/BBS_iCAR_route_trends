// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level abundance
// no spatial prior on the trends - for comparison with the spatial BYM trend model
// slope_iCAR_route2.stan

//iCAR function for abundance only
// functions {
//   real icar_normal_lpdf(vector bb, int nroutes, int[] node1, int[] node2) {
//     return -0.5 * dot_self(bb[node1] - bb[node2])
//       + normal_lpdf(sum(bb) | 0, 0.001 * nroutes); //soft sum to zero constraint on bb
//  }
// }


data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;
  int<lower=1> nstrata;
 
  int<lower=0> count[ncounts];              // count observations
  int<lower=1> year[ncounts]; // year index
  int<lower=1> route[ncounts]; // route index
  int<lower=0> firstyr[ncounts]; // first year index
  int<lower=1> observer[ncounts];              // observer indicators
  int<lower=1> strat[ncounts];              // observer indicators

  int<lower=1> fixedyear; // centering value for years
 
 // // spatial neighbourhood information
 //  int<lower=1> N_edges;
 //  int<lower=1, upper=nroutes> node1[N_edges];  // node1[i] adjacent to node2[i]
 //  int<lower=1, upper=nroutes> node2[N_edges];  // and node1[i] < node2[i]



}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion

  vector[nstrata] beta_raw_rand;
  real BETA; 

  vector[nroutes] alpha_raw;
  real STRATA; 

  vector[nstrata] strata_raw;
  
  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
 real<lower=0> sdbeta_rand;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts

  
}


model {


  vector[ncounts] E;           // log_scale additive likelihood
   vector[nstrata] beta_rand;
 vector[nstrata] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;
  
  vector[nstrata] strata;

// covariate effect on intercepts and slopes
   beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_rand + BETA;
   strata = strata_raw + STRATA;
   
   alpha = (sdalpha*alpha_raw);
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[strat[i]] * (year[i]-fixedyear) + strata[strat[i]] + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  
  beta_raw_rand ~ normal(0,1);//observer effects
  sum(beta_raw_rand) ~ normal(0,0.001*nroutes);

  strata_raw ~ normal(0,1);//observer effects
  sum(strata_raw) ~ normal(0,0.001*nroutes);

  
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1); //~ student_t(4,0,1); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
   
  
  alpha_raw ~ normal(0,1);//icar_normal(nroutes, node1, node2);
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  
  sdalpha ~ normal(0,1); //prior on sd of intercept variation
  

  sdbeta_rand  ~ gamma(2,10);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation



}

 generated quantities {

     vector[ncounts] log_lik;
     
       vector[ncounts] E;           // log_scale additive likelihood
   vector[nstrata] beta_rand;
  vector[nstrata] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;

  vector[nstrata] strata;
  
  matrix[nstrata,nyears] n;

// intercepts and slopes
  beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_rand + BETA;
      strata = strata_raw + STRATA;
 
    alpha = (sdalpha*alpha_raw);
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[strat[i]] * (year[i]-fixedyear) + strata[strat[i]] + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  
  
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  
  for(s in 1:nstrata){
    
    for(y in 1:nyears){
      n[s,y] = beta[s] * (y-fixedyear) + strata[s] + 0.5*sdalpha*sdalpha + + 0.5*sdobs*sdobs + 0.5*sdnoise*sdnoise;
    }
  }

 }

