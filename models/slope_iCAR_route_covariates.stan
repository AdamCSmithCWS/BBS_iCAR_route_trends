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
   int<lower=1> n_c_alpha; //number of intercept covariates
  int<lower=1> n_c_beta; //number of slope covariates
 
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

// Covariates
  //route-level covariates on slopes
  matrix[nroutes,n_c_beta] beta_covs;
  //route-level quadratic covariates on slopes
//  matrix[nroutes,n_c_beta] beta_covs2;

   //route-level covariates on intercepts
   matrix[nroutes,n_c_alpha] alpha_covs;
   //route-level quadratic covariates on intercepts
 //  matrix[nroutes,n_c_alpha] alpha_covs2;

}

parameters {

   vector[n_c_alpha] c_alpha; //covariate parameters on the intercepts
  // vector[n_c_alpha] c_alpha2; //covariate quadratic parameters on the intercepts

  vector[n_c_beta] c_beta;  // covariate parameters on the slopes
 // vector[n_c_beta] c_beta2;  // covariate quadratic parameters on the slopes

 
  vector[ncounts] noise_raw;             // over-dispersion

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


model {


  vector[ncounts] E;           // log_scale additive likelihood
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;
  vector[nroutes] sum_beta;
  vector[nroutes] sum_alpha;


//weakly informative normal priors assuming standardized predictors
    c_beta ~ normal(0,0.1);
//     c_beta2 ~ normal(0,0.05);
     

    sum_beta = beta_covs * c_beta;// + beta_covs2 * c_beta2;  //summary of the covariate effects on trends
  
     c_alpha ~ normal(0,0.1);
//     c_alpha2 ~ normal(0,0.05); //not sure why these polynomial terms aren't working...
     

    sum_alpha = alpha_covs * c_alpha;// + alpha_covs2 * c_alpha2;  //summary of the covariate effects on trends
 

  
  sdnoise ~ gamma(2,3);//boundary avoiding prior on sd
  // sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ student_t(4,0,1); //heavy tailed extra Poisson variance
  sum(noise_raw) ~ normal(0,0.0001*ncounts);
  
  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of slope variation

  beta_raw ~ icar_normal_lpdf(nroutes, node1, node2);
  alpha_raw ~ icar_normal_lpdf(nroutes, node1, node2);


// spatial and covariate effects on intercepts and slopes
   beta = (sdbeta*beta_raw) + sum_beta + BETA;
   alpha = (sdalpha*alpha_raw) + sum_alpha + ALPHA;
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;
 


  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
   count ~ poisson_log(E); //vectorized count likelihood with log-transformation
 

}

 generated quantities {

     vector[ncounts] log_lik;
     
       vector[ncounts] E;           // log_scale additive likelihood
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;
  vector[nroutes] sum_beta;
  vector[nroutes] sum_alpha;

    sum_beta = beta_covs * c_beta ;//+ beta_covs2 * c_beta2;  //summary of the covariate effects on trends
 sum_alpha = alpha_covs * c_alpha ;//+ alpha_covs2 * c_alpha2;  //summary of the covariate effects on trends
  
// covariate effect on intercepts and slopes
   beta = (sdbeta*beta_raw) + sum_beta + BETA;
   alpha = (sdalpha*alpha_raw) + sum_alpha + ALPHA;
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  
  
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  

 }

