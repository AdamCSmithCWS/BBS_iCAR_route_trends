// This is a Stan implementation of a route-level slope model
// plus, it has an explicitly spatial prior structure on the 
// random effect, stratum-level trends
// and no random year-effects - slope only


//iCAR function
functions {
 //  real icar_normal_lpdf(vector bb, int nroutes, int[] node1, int[] node2) {
 //    return -0.5 * dot_self(bb[node1] - bb[node2])
 //      + normal_lpdf(sum(bb) | 0, 0.001 * nroutes); //soft sum to zero constraint on bb
 // }
 // 
 
 /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
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
 
  matrix<lower = 0, upper = 1>[nroutes, nroutes] W; // adjacency matrix
  int W_n;                // number of adjacent region pairs = sum(W)/2
  
  

}


transformed data {
  // transformed information for spatial components
    int W_sparse[W_n, 2];   // adjacency pairs
  vector[nroutes] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[nroutes] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(nroutes - 1)) {
      for (j in (i + 1):nroutes) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:nroutes) D_sparse[i] = sum(W[i]);
  {
    vector[nroutes] invsqrtD;  
    for (i in 1:nroutes) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}



parameters {
  vector[ncounts] noise_raw;             // over-dispersion

  vector[nroutes] beta_space;
  vector[nroutes] beta_raw_rand;
  real BETA; 

  vector[nroutes] alpha_space;
  //real ALPHA; 

  real eta; //first-year intercept
  
  vector[nobservers] obs_raw; //observer effects

  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta_space;    // sd of slopes 
 real<lower=0> sdbeta_rand;    // sd of slopes 
  real<lower=0> sdalpha_space;    // sd of intercepts


  //real<lower = 0> tau_alpha; //
  real<lower = 0, upper = 1> a_alpha; //spatial covariance for abundance
  //real<lower = 0> tau_beta;
  real<lower = 0, upper = 1> a_beta; //spatial covariance for slopes
  
  
}

model {


  vector[ncounts] E;           // log_scale additive likelihood
   vector[nroutes] beta_rand;
  //vector[nroutes] beta_space;
 vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;

// covariate effect on intercepts and slopes
   //beta_space = (sdbeta_space*beta_raw_space);
   beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_space + beta_rand + BETA;
   alpha = alpha_space;// + ALPHA;
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  
  beta_raw_rand ~ normal(0,1);//random non-spatial variation in slope
  sum(beta_raw_rand) ~ normal(0,0.001*nroutes);

  
  sdnoise ~ normal(0,0.5); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1); //~ student_t(4,0,1); //normal tailed extra Poisson log-normal variance
  
  sdobs ~ std_normal(); //prior on sd of gam hyperparameters
 
  obs_raw ~ normal(0,1);//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  //ALPHA ~ normal(0,1);// prior on fixed effect mean intercept
  eta ~ normal(0,1);// prior on first-year observer effect
  a_beta ~ uniform(0,1);//
  a_alpha ~ uniform(0,1);//

  //spatial CAR intercepts and slopes by strata
  sdalpha_space ~ gamma(2,2);//~ normal(0,1); //prior on sd of intercept variation
  sdbeta_space ~ gamma(2,20);//~ normal(0,0.05); //boundary avoiding prior on sd of slope spatial variation w mean = 0.1 and 99% < 0.33
  sdbeta_rand  ~ gamma(2,20);//~ normal(0,0.05); //boundary avoiding prior on sd of slope random variation

  beta_space ~ sparse_car(sdbeta_space, a_beta, W_sparse, D_sparse, lambda, nroutes, W_n);
  alpha_space ~ sparse_car(sdalpha_space, a_alpha, W_sparse, D_sparse, lambda, nroutes, W_n);


// vector phi, real tau, real alpha, 
//     int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n
//     
}

 generated quantities {

     vector[ncounts] log_lik;
     
       vector[ncounts] E;           // log_scale additive likelihood
   vector[nroutes] beta_rand;
  //vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  vector[ncounts] noise;

// intercepts and slopes
   //beta_space = (sdbeta_space*beta_raw_space);
   beta_rand = (sdbeta_rand*beta_raw_rand);
   
   beta = beta_space + beta_rand + BETA;
    alpha = alpha_space;// + ALPHA;
   noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i] + noise[i];
  }
  
  
  
  for(i in 1:ncounts){
  log_lik[i] = poisson_log_lpmf(count[i] | E[i]);
  }
  

 }
