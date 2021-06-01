### a function to calculate the basis function(s) for a GAM component of a stan model

gam_basis <- function(orig.preds = dts[,"yr"],
                           nknots = 6,#number of internal knots
                           predpoints = NULL,
                           npredpoints = 100,
                           sm_name = ""){
  
  if(any(is.na(orig.preds) == T)){
    stop("This GAM formulation cannot handle missing values in the predictor")
  }
  
  
  require(mgcv)
    
    dat = data.frame(x = orig.preds,
                     y = rnorm(length(orig.preds),0,0.1))
    if(is.null(predpoints)){
    predpoints = seq(min(orig.preds),max(orig.preds),length.out = npredpoints)
    }else{
      npredpoints <- length(predpoints)
    }
    
    dat_pred = data.frame(x = predpoints,
                          y = rnorm(length(predpoints),0,0.1))
    

    M = smoothCon(s(x,k = nknots+1, bs = "tp"),data = dat,
                       absorb.cons=TRUE,#this drops the constant
                       diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 
    
    
    gamx.basis = M[[1]]$X
    
    
    M_pred = smoothCon(s(x,k = nknots+1, bs = "tp"),data = dat_pred,
                       absorb.cons=TRUE,#this drops the constant
                       diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 
    
    gamx.basispred = M_pred[[1]]$X
    
    
    
    
   outlist <- list(gamx.basis = gamx.basis,
                  gamx.basispred = gamx.basispred,
                  orig.preds = orig.preds,
                  predpoints = predpoints,
                  nknots = nknots,
                  npredpoints = npredpoints,
                  M = M,
                  M_pred = M_pred)
  names(outlist) <- c(paste0(sm_name,"_basis"),
                      paste0(sm_name,"_basispred"),
                      "original_predictor_values",
                      paste0(sm_name,"_visualized_predictor_values"),
                      paste0("nknots_",sm_name),
                      paste0("npredpoints_",sm_name),
                      paste0(sm_name,"_smoothCon"),
                      paste0(sm_name,"_smoothCon_pred"))
  
  
  return(outlist)
  
  
  
}
