

### function to extract the dimension values from an bayesian model fit
### works within the gather_samples function
dim_ext <- function(dim = 1,
                     var = "",
                     cl = "Parameter",
                     dat = NULL){
  ##3 function to extract the indicator values from cmdstanr output
  require(stringr)
  
  pat = paste0("(?<=",var,"\\[")
  
  if(dim > 1){
    for(j in 1:(dim-1)){
      
      pat2 = paste0(pat,")[:digit:]+")
      cl2 = str_extract(unlist(dat[,cl]),pattern = pat2)
      
      d = max(nchar(cl2))
      
      pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
    }
  }
  
  
  pat = paste0(pat,")[:digit:]+")
  dds = as.integer(str_extract(unlist(dat[,cl]),pattern = pat))
  return(dds)
  
}


### function to generate the same tidy output as gather-draws in tidybayes package
## dims should be a character vector defining the dimensions of the parameter
## e.g., parm = "nsmooth", dims = c("stratum","year"),
## function works with cmdstanr output by default and rstan fits
## with is_rstan == TRUE
posterior_samples <- function(fit = cmdstanfit,
                        parm = "nsmooth",
                        dims = NULL,
                        is_rstan = FALSE,
                        is_mcmc = FALSE){
  require(posterior)
  require(tidyverse)
  
  if(length(dims) > 0){
    parm_ex <- paste0(parm,"[")
  }else{
    parm_ex <- parm
  }
  if(class(fit)[1] == "stanfit"){is_rstan <- TRUE}
  
  if(class(fit)[1] == "mcmc"){is_mcmc <- TRUE}
  
  if(is_rstan | is_mcmc){
  samples <- as_draws_df(as.array(fit)) %>% 
    dplyr::select(starts_with(parm_ex,ignore.case = FALSE),
           .chain,
           .iteration,
           .draw)#,pars = c(parm)))
  }else{

    samples <- as_draws_df(fit$draws(variables = c(parm)))
    
  }
  

  
  plong <- suppressWarnings(samples %>% pivot_longer(
    cols = starts_with(parm_ex,ignore.case = FALSE),
    names_to = c(".variable"),
    values_to = ".value",
    values_drop_na = TRUE
  )) 
  
  for(dn in 1:length(dims)){
    dd = dims[dn]
    plong[,dd] = dim_ext(dim = dn,
                          var = parm,
                          cl = ".variable",
                          dat = plong)
    
  }
  
  plong <- plong %>% mutate(.variable = parm)
  return(plong)
  
}





posterior_sums <- function(samples = n_samples,
                      quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                      ci = 0.95,
                      dims = NULL){
  
  cis = c((1-ci)/2,1-((1-ci)/2))
 
   if(!is.null(dims)){
     for(i in 1:length(dims)){
       samples <- samples %>% 
         rename_with(~gsub(dims[i],paste0("d",i),.x,fixed = TRUE))
     }
    if(length(dims) == 1){
     sums = samples %>% 
       group_by(d1) %>% 
       summarise(mean = mean(.value),
                 median = median(.value),
                 sd = sd(.value),
                 lci = quantile(.value,cis[1]),
                 uci = quantile(.value,cis[2])) %>% 
       rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))
     
     if(!is.null(quantiles)){
       qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
       for(i in 1:length(quantiles)){
         qq = quantiles[i]
         qn = qs[i]
         
         sumt = samples %>% 
           group_by(d1) %>% 
           summarise(tt = as.numeric(quantile(.value,qq))) %>% 
           rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
           rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) 
           
         sums = left_join(sums,sumt,by = dims)
       }
       
     }
     
    }
     if(length(dims) == 2){
       sums = samples %>% 
         group_by(d1,d2) %>% 
         summarise(mean = mean(.value),
                   median = median(.value),
                   sd = sd(.value),
                   lci = quantile(.value,cis[1]),
                   uci = quantile(.value,cis[2])) %>% 
         rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
         rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))

       
       if(!is.null(quantiles)){
         qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
         for(i in 1:length(quantiles)){
           qq = quantiles[i]
           qn = qs[i]
           
           sumt = samples %>% 
             group_by(d1,d2) %>% 
             summarise(tt = as.numeric(quantile(.value,qq))) %>% 
             rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
             rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
             rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))
           
           
           sums = left_join(sums,sumt,by = dims)
         }
         
       }
       
       
       }
     
     if(length(dims) == 3){
       sums = samples %>% 
         group_by(d1,d2,d3) %>% 
         summarise(mean = mean(.value),
                   median = median(.value),
                   sd = sd(.value),
                   lci = quantile(.value,cis[1]),
                   uci = quantile(.value,cis[2])) %>% 
         rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
         rename_with(~gsub("d2",dims[2],.x,fixed = TRUE)) %>% 
         rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))
       
       
       if(!is.null(quantiles)){
         qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
         for(i in 1:length(quantiles)){
           qq = quantiles[i]
           qn = qs[i]
           
           sumt = samples %>% 
             group_by(d1,d2,d3) %>% 
             summarise(tt = as.numeric(quantile(.value,qq))) %>% 
             rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
             rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
             rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))  %>% 
             rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))
           
           
           sums = left_join(sums,sumt,by = dims)
         }
         
       }
       
       
     }
     

     if(length(dims) == 4){
       sums = samples %>% 
         group_by(d1,d2,d3,d4) %>% 
         summarise(mean = mean(.value),
                   median = median(.value),
                   sd = sd(.value),
                   lci = quantile(.value,cis[1]),
                   uci = quantile(.value,cis[2])) %>% 
         rename_with(~gsub("d1",dims[1],.x,fixed = TRUE)) %>% 
         rename_with(~gsub("d2",dims[2],.x,fixed = TRUE)) %>% 
         rename_with(~gsub("d3",dims[3],.x,fixed = TRUE)) %>% 
         rename_with(~gsub("d4",dims[4],.x,fixed = TRUE))
       
       
       if(!is.null(quantiles)){
         qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
         for(i in 1:length(quantiles)){
           qq = quantiles[i]
           qn = qs[i]
           
           sumt = samples %>% 
             group_by(d1,d2,d3,d4) %>% 
             summarise(tt = as.numeric(quantile(.value,qq))) %>% 
             rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE)) %>% 
             rename_with(~gsub("d1",dims[1],.x,fixed = TRUE))  %>% 
             rename_with(~gsub("d2",dims[2],.x,fixed = TRUE))  %>% 
             rename_with(~gsub("d3",dims[3],.x,fixed = TRUE))  %>% 
             rename_with(~gsub("d4",dims[4],.x,fixed = TRUE))
           
           
           sums = left_join(sums,sumt,by = dims)
         }
         
       }
       
       
     }
     
     if(length(dims) > 4){stop("ERROR this function cannot handle more than 4 dimensions, but it could be easily modified")}
    
  }else{
  
 
  sums = samples %>% 
    summarise(mean = mean(.value),
              median = median(.value),
              sd = sd(.value),
              lci = quantile(.value,cis[1]),
              uci = quantile(.value,cis[2]))
  
  if(!is.null(quantiles)){
    qs = paste0("Q_",gsub(quantiles,pattern = "0.",replacement = "",fix = TRUE))
    for(i in 1:length(quantiles)){
      qq = quantiles[i]
      qn = qs[i]
      
      sumt = samples %>% 
        summarise(tt = as.numeric(quantile(.value,qq))) %>% 
        rename_with(~gsub(pattern = "tt",replacement = qn,.x,fixed = TRUE))
      sums = bind_cols(sums,sumt)
    }
    
  }
  
  }
  return(sums)
  
}



