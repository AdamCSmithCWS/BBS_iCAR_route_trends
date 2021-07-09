// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only



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
  

  
  // values for predicting next years data
  int<lower=1> ncounts_pred;
  int<lower=0> count_pred[ncounts_pred];              // count observations
  int<lower=1> route_pred[ncounts_pred]; // route index
  int<lower=0> firstyr_pred[ncounts_pred]; // first year index
  int<lower=1> observer_pred[ncounts_pred]; 
  
}

parameters {
  vector[ncounts] noise_raw;             // over-dispersion
  
  vector[nroutes] beta_raw_rand;
  real BETA; 
  
  vector[nroutes] alpha_raw;
  real ALPHA; 
  
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
  vector[nroutes] beta_rand;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;
  
  // covariate effect on intercepts and slopes
  beta_rand = (sdbeta_rand*beta_raw_rand);
  
  beta = beta_rand + BETA;
  alpha = (sdalpha*alpha_raw) + ALPHA;
  noise = sdnoise*noise_raw;
  obs = sdobs*obs_raw;
  
  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  
  beta_raw_rand ~ normal(0,1);//observer effects
  sum(beta_raw_rand) ~ normal(0,0.001*nroutes);
  
  
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1); //~ student_t(4,0,1); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
  
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);
  
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,1); //prior on sd of intercept variation
  sdbeta_rand  ~ gamma(2,20);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation
  
  alpha_raw ~ normal(0,1);//icar_normal(nroutes, node1, node2);
  sum(alpha_raw) ~ normal(0,0.001*nroutes);
  
  
}

generated quantities {
  
  vector[ncounts_pred] log_lik;
  
  vector[ncounts_pred] E_pred;           // log_scale additive likelihood
  vector[nroutes] beta_rand;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  
  // intercepts and slopes
  beta_rand = (sdbeta_rand*beta_raw_rand);
  
  beta = beta_rand + BETA;
  alpha = (sdalpha*alpha_raw) + ALPHA;
  obs = sdobs*obs_raw;
  
  // Predictions for next year's data

  for(i in 1:ncounts_pred){
    real obs_tmp;
    real noise_tmp = sdnoise*normal_rng(0,1);
    if(firstyr_pred[i] == 0) 
    obs_tmp = obs[observer_pred[i]];
    else
    obs_tmp = sdobs*normal_rng(0,1);
    
    E_pred[i] =  beta[route_pred[i]] * ((nyears+1)-fixedyear) + alpha[route_pred[i]] + obs_tmp + eta*firstyr_pred[i] + noise_tmp;
  }
  
  
  
  for(i in 1:ncounts_pred){
  log_lik[i] = poisson_log_lpmf(count_pred[i] | E_pred[i]);
  }
  

 }

