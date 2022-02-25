# initial values for different models

init_def_BYM <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                             alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                             ALPHA = 0,
                             BETA = 0,
                             eta = 0,
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             sdnoise = 0.2,
                             sdobs = 0.1,
                             sdbeta_space = runif(1,0.01,0.1),
                             sdbeta_rand = runif(1,0.01,0.1),
                             beta_raw_space = rnorm(stan_data$nroutes,0,0.01),
                             beta_raw_rand = rnorm(stan_data$nroutes,0,0.01))}

init_def_iCAR <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                                 alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                                 ALPHA = 0,
                                 BETA = 0,
                                 eta = 0,
                                 obs_raw = rnorm(stan_data$nobservers,0,0.1),
                                 sdnoise = 0.2,
                                 sdobs = 0.1,
                                 sdbeta_space = runif(1,0.01,0.1),
                                 beta_raw_space = rnorm(stan_data$nroutes,0,0.01))}

init_def_Non_spatial <- function(){ list(noise_raw = rnorm(stan_data$ncounts,0,0.1),
                                 alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                                 ALPHA = 0,
                                 BETA = 0,
                                 eta = 0,
                                 obs_raw = rnorm(stan_data$nobservers,0,0.1),
                                 sdnoise = 0.2,
                                 sdobs = 0.1,
                                 sdbeta_rand = runif(1,0.01,0.1),
                                 beta_raw_rand = rnorm(stan_data$nroutes,0,0.01))}

init_func_list <- list(BYM = init_def_BYM,
                       iCAR = init_def_iCAR,
                       Non_spatial = init_def_Non_spatial)

