// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only


//iCAR function
 functions {
   real icar_normal_lpdf(vector bb, int ns, array[] int n1, array[] int n2) {
     return -0.5 * dot_self(bb[n1] - bb[n2])
       + normal_lpdf(sum(bb) | 0, 0.001 * ns); //soft sum to zero constraint on bb
  }
 }

data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> nobservers;
 
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]


}

parameters {

  vector[nroutes] beta_raw_space;
  vector[nroutes] beta_raw_rand;
  real BETA; 

  vector[nroutes] alpha_raw;
  real ALPHA; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta_space;    // sd of slopes 
 real<lower=0> sdbeta_rand;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts

  
}


model {


  vector[ncounts] E;           // log_scale additive likelihood
   vector[nroutes] beta_rand;
  vector[nroutes] beta_space;
 vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  real phi;

// covariate effect on intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_space + beta_rand + BETA;
   alpha = (sdalpha*alpha_raw) + ALPHA;
 //  noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i];
  }
  
  
  beta_raw_rand ~ normal(0,1);//random slope effects
  sum(beta_raw_rand) ~ normal(0,0.001*nroutes);

  
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
//  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  //spatial iCAR intercepts and slopes by strata
  sdalpha ~ normal(0,2); //prior on sd of intercept variation
  sdbeta_space ~ gamma(2,50);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.04 and 99% < 0.13
  sdbeta_rand  ~ gamma(2,50);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation

  beta_raw_space ~ icar_normal(nroutes, node1, node2);
  alpha_raw ~ icar_normal(nroutes, node1, node2);


}

 generated quantities {

   vector[nroutes] beta_rand;
  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;

// intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_space + beta_rand + BETA;
    alpha = (sdalpha*alpha_raw) + ALPHA;
  

 }

