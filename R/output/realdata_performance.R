####################################################################################
###
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview:
###
####################################################################################

rm(list = ls())

#########################################################################

set.seed(10112021)

library(dplyr)
library(GpGp)
library(correlationVecchia)

#########################################################################

# b2 = S-E-MM + J-E-NN
# b3 = S-E-MM + S-E-NN
# b5 = T-ord + T-NN
# b7 = T-ord + J-C-NN
# b8 = T-ord + S-C-NN
# cc = C-MM + C-NN

load("DATA/joint_CRCM_NCEP_10152021_estimates.RData")

round(fit.m50[[1]]$betahat, 5) ; round(fit.m50[[1]]$covparms, 5) # S-E-MM + J-E-NN
round(fit.m50[[2]]$betahat, 5) ; round(fit.m50[[2]]$covparms, 5) # S-E-MM + S-E-NN
round(fit.m50[[3]]$betahat, 5) ; round(fit.m50[[3]]$covparms, 5) # T-ord + T-NN
round(fit.m50[[4]]$betahat, 5) ; round(fit.m50[[4]]$covparms, 5) # T-ord + J-C-NN
round(fit.m50[[5]]$betahat, 5) ; round(fit.m50[[5]]$covparms, 5) # T-ord + S-C-NN
round(fit.m50[[6]]$betahat, 5) ; round(fit.m50[[6]]$covparms, 5) # C-MM + C-NN

ms          <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
fit.all     <- list(m5 = fit.m5, m10 = fit.m10, m15 = fit.m15, m20 = fit.m20, m25 = fit.m25, m30 = fit.m30, m35 = fit.m35, m40 = fit.m40, m45 = fit.m45, m50 = fit.m50)

rm(fit.m5, fit.m10, fit.m15, fit.m20, fit.m25, fit.m30, fit.m35, fit.m40, fit.m45, fit.m50)

#########################################################################

table.runtime <- matrix(NA, nrow = length(fit.all), ncol = 6)

for(i in 1:length(fit.all)) {
  for(j in 1:6) {
    
    table.runtime[i, j] <- fit.all[[i]][[6 + j]][3]
  }
}

table.runtime <- as.data.frame(table.runtime)

rownames(table.runtime) <- paste0("m = ", ms)
colnames(table.runtime) <- c("S-E-MM + J-E-NN", "S-E-MM + S-E-NN", "T-ord + T-NN", "T-ord + J-C-NN", "T-ord + S-C-NN", "C-MM + C-NN")

# table.runtime[, c(2, 1, 3, 5, 4, 6), drop = FALSE]

#########################################################################

colnames(df.joint)

locs        <- train.joint[, c("x1", "x2", "t", "d")] %>% as.matrix()
z           <- train.joint$z
X.train     <- cbind(1, train.joint$d) %>% as.matrix()

locs.pred   <- test.joint[, c("x1", "x2", "t", "d")] %>% as.matrix()
z.pred      <- test.joint$z
X.pred      <- cbind(1, test.joint$d) %>% as.matrix()

ns_obs      <- train.joint %>% pull(d) %>% table() %>% as.numeric()
ns_pred     <- test.joint %>% pull(d) %>% table() %>% as.numeric()

rm(train.joint, test.joint)

#########################################################################

predictions_scaled_modified_for_b3 <- function(fit, locs_pred, m, joint, nsims, predvar, X_pred, scale, tol)
{
  
  y_obs = fit$y
  locs_obs = fit$locs
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf
  
  # ## add nugget for numerical stability
  # if(covparms[length(covparms)]==0)
  #   covparms[length(covparms)]=covparms[1]*1e-12
  
  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }
  
  # specify how to scale input dimensions
  if(scale=='parms'){ scales=1/covparms[1+(1:ncol(locs_obs))]
  } else if(scale=='ranges'){ scales=1/apply(locs_obs,2,function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=',scale))
  
  
  ###
  if(joint){  # joint predictions
    
    # get orderings
    temp=order_maxmin_pred(locs_obs,locs_pred)
    ord1=temp$ord
    ord2=temp$ord_pred
    
    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta
    
    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    
    # get nearest neighbor array
    sm = if (n_pred<1e5) 2 else 1.5
    NNarray_all <- find_ordered_nn_pred(locs_all,m,fix.first=n_obs,searchmult=sm)
    NNarray_pred=NNarray_all[-(1:n_obs),-1]
    
    means=numeric(length=n_pred)
    if(nsims>0) samples=array(dim=c(n_pred,nsims))
    
    # make predictions sequentially
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=sort(NNarray_pred[i,])
      NN_obs=NN[NN<=n_obs]
      NN_pred=NN[NN>n_obs]-n_obs
      
      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(locs_all[c(NN,i+n_obs),]))
      cl=t(chol(K + tol * diag(nrow(K))))
      
      # prediction
      y.NN=y[NN_obs]
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,means[NN_pred]))
      if(nsims>0){ # conditional simulation
        pred.var=cl[m+1,m+1]^2*vcf
        for(s in 1:nsims){
          pm=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,samples[NN_pred,s]))
          samples[i,s]=stats::rnorm(1,pm,sqrt(pred.var))
        }
      }
      
    }
    
    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta)
    if(nsims==0){
      preds=means
    } else {
      samples[ord2,] = samples + c(Xord_pred %*% beta)
      preds=list(means=means,samples=samples)
    }
    
  } else {  # separate predictions
    
    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')
    
    y  = y_obs - X_obs %*% beta
    
    # find the NNs
    m=min(m,nrow(locs_obs))
    NNarray=FNN::get.knnx(locs_obs,locs_pred,m)$nn.index
    
    
    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=NNarray[i,]
      
      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(rbind(locs_obs[NN,],locs_pred[i,])))
      cl=t(chol(K + tol * diag(nrow(K))))
      
      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf
      
    }
    means=means+c(X_pred %*% beta)
    
    if(predvar==FALSE){
      preds=means
    } else {
      preds=list(means=means,vars=vars)
    }
    
  }
  
  return(preds)
}

predictions_scaled_modified_for_b8 <- function(fit, locs_pred, m, joint, nsims, predvar, X_pred, scale, tol)
{
  
  y_obs = fit$y
  locs_obs = fit$locs
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  if(is.null(fit$vcf)) vcf=1 else vcf=fit$vcf
  
  # ## add nugget for numerical stability
  # if(covparms[length(covparms)]==0)
  #   covparms[length(covparms)]=covparms[1]*1e-12
  
  # specify trend if missing
  if(missing(X_pred)){
    if(fit$trend=='zero'){
      X_pred=as.matrix(rep(0,n_pred))
    } else if(fit$trend=='intercept'){
      X_pred=as.matrix(rep(1,n_pred))
    } else if(fit$trend=='linear'){
      X_pred=cbind(rep(1,n_pred),locs_pred)
    } else stop('X_pred must be specified')
  }
  
  # specify how to scale input dimensions
  if(scale=='parms'){ scales=1/covparms[1+(1:ncol(locs_obs))]
  } else if(scale=='ranges'){ scales=1/apply(locs_obs,2,function(x) diff(range(x)))
  } else stop(paste0('invalid argument scale=',scale))
  
  
  ###
  if(joint){  # joint predictions
    
    # get orderings
    temp=list(ord = order_time(locs_obs), ord_pred = order_time(locs_pred))
    ord1=temp$ord
    ord2=temp$ord_pred
    
    # reorder stuff
    X_obs <- as.matrix(X_obs)
    X_pred <- as.matrix(X_pred)
    Xord_obs  <- X_obs[ord1,,drop=FALSE]
    Xord_pred <- X_pred[ord2,,drop=FALSE]
    y  <- y_obs[ord1] - Xord_obs %*% beta
    
    # put all coordinates together
    locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
    
    # get nearest neighbor array
    sm = if (n_pred<1e5) 2 else 1.5
    NNarray_all <- find_ordered_nn_pred(locs_all,m,fix.first=n_obs,searchmult=sm)
    NNarray_pred=NNarray_all[-(1:n_obs),-1]
    
    means=numeric(length=n_pred)
    if(nsims>0) samples=array(dim=c(n_pred,nsims))
    
    # make predictions sequentially
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=sort(NNarray_pred[i,])
      NN_obs=NN[NN<=n_obs]
      NN_pred=NN[NN>n_obs]-n_obs
      
      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(locs_all[c(NN,i+n_obs),]))
      cl=t(chol(K + tol * diag(nrow(K))))
      
      # prediction
      y.NN=y[NN_obs]
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,means[NN_pred]))
      if(nsims>0){ # conditional simulation
        pred.var=cl[m+1,m+1]^2*vcf
        for(s in 1:nsims){
          pm=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],c(y.NN,samples[NN_pred,s]))
          samples[i,s]=stats::rnorm(1,pm,sqrt(pred.var))
        }
      }
      
    }
    
    # add (prior) mean and return to original ordering
    means[ord2] = means + c(Xord_pred %*% beta)
    if(nsims==0){
      preds=means
    } else {
      samples[ord2,] = samples + c(Xord_pred %*% beta)
      preds=list(means=means,samples=samples)
    }
    
  } else {  # separate predictions
    
    if(nsims>0) stop('cannot produce joint samples when joint=FALSE')
    
    y  = y_obs - X_obs %*% beta
    
    # find the NNs
    m=min(m,nrow(locs_obs))
    NNarray=FNN::get.knnx(t(t(locs_obs)*scales),
                          t(t(locs_pred)*scales),m)$nn.index
    
    
    means=vars=numeric(length=n_pred)
    for(i in 1:n_pred){
      
      # NN conditioning sets
      NN=NNarray[i,]
      
      # (co-)variances (modified)
      K=get(covfun_name)(covparms, as.matrix(rbind(locs_obs[NN,],locs_pred[i,])))
      cl=t(chol(K + tol * diag(nrow(K))))
      
      # prediction
      means[i]=cl[m+1,1:m]%*%forwardsolve(cl[1:m,1:m],y[NN])
      vars[i]=cl[m+1,m+1]^2*vcf
      
    }
    means=means+c(X_pred %*% beta)
    
    if(predvar==FALSE){
      preds=means
    } else {
      preds=list(means=means,vars=vars)
    }
    
  }
  
  return(preds)
}

#########################################################################

spliting.fit <- function(fit)
{
  fit$covparms <- fit$covparms[-5]
  fit$logparms <- fit$logparms[-5]
  fit$grad <- fit$grad[-5]
  fit$info <- fit$info[-5, -5]
  
  # fit$loglik <- fit$loglik
  # fit$no_decrease <- fit$no_decrease
  # fit$conv <- fit$conv
  # fit$covfun_name <- fit$covfun_name
  # fit$trend <- fit$trend
  # fit$vcf <- fit$vcf
  
  ##
  
  idx.d0 <- which(fit$locs[,4] == 0)
  idx.d1 <- which(fit$locs[,4] == 1)
  
  fit.d0 <- fit
  fit.d1 <- fit
  
  ##
  
  fit.d0$betahat <- fit$betahat[1]
  fit.d0$sebeta <- NA
  fit.d0$tbeta <- NA
  fit.d0$betacov <- NA
  fit.d0$y <- fit.d0$y[idx.d0]
  fit.d0$locs <- fit.d0$locs[idx.d0, -4, drop = FALSE]
  fit.d0$X <- fit.d0$X[idx.d0, 1, drop = FALSE] 
  
  ##
  
  fit.d1$betahat <- fit$betahat[1] + fit$betahat[2]
  fit.d1$sebeta <- NA
  fit.d1$tbeta <- NA
  fit.d1$betacov <- NA
  fit.d1$y <- fit.d1$y[idx.d1]
  fit.d1$locs <- fit.d1$locs[idx.d1, -4, drop = FALSE]
  fit.d1$X <- fit.d1$X[idx.d1, 1, drop = FALSE] 
  
  ##
  
  return( list(d0 = fit.d0, d1 = fit.d1) )
}

predictions_split <- function(approx, fit, locs_pred, m, joint, nsims, predvar, X_pred, scale, tol)
{
  ##
  
  fit <- spliting.fit(fit)
  idx.d0 <- which(locs_pred[, 4] == 0)
  idx.d1 <- which(locs_pred[, 4] == 1)
  
  locs_pred.d0 <- locs_pred[idx.d0, -4, drop = FALSE]
  locs_pred.d1 <- locs_pred[idx.d1, -4, drop = FALSE]
  
  X_pred.d0 <- X_pred[idx.d0, 1, drop = FALSE]
  X_pred.d1 <- X_pred[idx.d1, 1, drop = FALSE]
  
  ##
  
  if(approx == 3) {
    
    pre.d0 <- predictions_scaled_modified_for_b3(fit = fit$d0, locs_pred = locs_pred.d0, m = m, joint = joint, nsims = nsims, predvar = predvar, X_pred = X_pred.d0, scale = scale, tol = tol)
    pre.d1 <- predictions_scaled_modified_for_b3(fit = fit$d1, locs_pred = locs_pred.d1, m = m, joint = joint, nsims = nsims, predvar = predvar, X_pred = X_pred.d1, scale = scale, tol = tol)
    
  } else if(approx == 8) {
    
    pre.d0 <- predictions_scaled_modified_for_b8(fit = fit$d0, locs_pred = locs_pred.d0, m = m, joint = joint, nsims = nsims, predvar = predvar, X_pred = X_pred.d0, scale = scale, tol = tol)
    pre.d1 <- predictions_scaled_modified_for_b8(fit = fit$d1, locs_pred = locs_pred.d1, m = m, joint = joint, nsims = nsims, predvar = predvar, X_pred = X_pred.d1, scale = scale, tol = tol)
    
  } else {
    
    stop("This function is only for b3 and b8. Check the argument approx. It must be 3 or 8.")
  }
  
  ##
  
  if(predvar == FALSE){
    
    means <- rep(NA, length(idx.d0) + length(idx.d1))
    means[idx.d0] <- pre.d0
    means[idx.d1] <- pre.d1
    
    preds <- means
    
  } else {
    
    means <- rep(NA, length(idx.d0) + length(idx.d1))
    means[idx.d0] <- pre.d0$means
    means[idx.d1] <- pre.d1$means
    
    vars <- rep(NA, length(idx.d0) + length(idx.d1))
    vars[idx.d0] <- pre.d0$vars
    vars[idx.d1] <- pre.d1$vars
    
    preds <- list(means = means, vars = vars)
  }
  
  ##
  
  return( preds )
}

#########################################################################

comparison <- function(m, joint = FALSE, idx = NULL)
{
  ##
  
  if(joint == TRUE) {
    
    predvar = FALSE 
    
  } else {
    
    predvar = TRUE
  }
  
  ##
  
  if(is.null(idx)) idx <- which(ms == m)

  fit.b2  <- fit.all[[idx]][[1]]
  fit.b3  <- fit.all[[idx]][[2]]
  fit.b5  <- fit.all[[idx]][[3]]
  fit.b7  <- fit.all[[idx]][[4]]
  fit.b8  <- fit.all[[idx]][[5]]
  fit.cc  <- fit.all[[idx]][[6]]
  
  ##
  
  # b2 = S-E-MM + J-E-NN
  # b3 = S-E-MM + S-E-NN
  # b5 = T-ord + T-NN
  # b7 = T-ord + J-C-NN
  # b8 = T-ord + S-C-NN
  # cc = C-MM + C-NN
  
  pre.b2  <- predictions_bs_mulv(approx = 2, ns_obs, ns_pred, fit = fit.b2, locs_pred = locs.pred, m = m, joint = joint, nsims = 0, predvar = predvar, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b3  <- predictions_split(  approx = 3,                  fit = fit.b3, locs_pred = locs.pred, m = m, joint = joint, nsims = 0, predvar = predvar, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b5  <- predictions_bs_mulv(approx = 5, ns_obs, ns_pred, fit = fit.b5, locs_pred = locs.pred, m = m, joint = joint, nsims = 0, predvar = predvar, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b7  <- predictions_bs_mulv(approx = 7, ns_obs, ns_pred, fit = fit.b7, locs_pred = locs.pred, m = m, joint = joint, nsims = 0, predvar = predvar, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.b8  <- predictions_split(  approx = 8,                  fit = fit.b8, locs_pred = locs.pred, m = m, joint = joint, nsims = 0, predvar = predvar, X_pred = X.pred, scale = 'parms', tol = 1e-8)
  pre.cc  <- predictions_cv_mulv(            ns_obs, ns_pred, fit = fit.cc, locs_pred = locs.pred, m = m, joint = joint, nsims = 0, predvar = predvar, X_pred = X.pred, scale = 'parms', tol = 1e-8)

  ##
  
  if(joint == TRUE) {
    
    mspe.b2 <- mean( (pre.b2 - z.pred)^2 )
    mspe.b3 <- mean( (pre.b3 - z.pred)^2 )
    mspe.b5 <- mean( (pre.b5 - z.pred)^2 )
    mspe.b7 <- mean( (pre.b7 - z.pred)^2 )
    mspe.b8 <- mean( (pre.b8 - z.pred)^2 )
    mspe.cc <- mean( (pre.cc - z.pred)^2 )
    
    return(list(mspe = c(mspe.b2, mspe.b3, mspe.b5, mspe.b7, mspe.b8, mspe.cc), mls = c(NA, NA, NA, NA, NA, NA)))
    
  } else if(joint == FALSE) {
   
    mspe.b2 <- mean( (pre.b2$means - z.pred)^2 )
    mspe.b3 <- mean( (pre.b3$means - z.pred)^2 )
    mspe.b5 <- mean( (pre.b5$means - z.pred)^2 )
    mspe.b7 <- mean( (pre.b7$means - z.pred)^2 )
    mspe.b8 <- mean( (pre.b8$means - z.pred)^2 )
    mspe.cc <- mean( (pre.cc$means - z.pred)^2 )
    
    mls.b2  <- mean(correlationVecchia::logscore(pre.b2$means, pre.b2$vars, z.pred))
    mls.b3  <- mean(correlationVecchia::logscore(pre.b3$means, pre.b3$vars, z.pred))
    mls.b5  <- mean(correlationVecchia::logscore(pre.b5$means, pre.b5$vars, z.pred))
    mls.b7  <- mean(correlationVecchia::logscore(pre.b7$means, pre.b7$vars, z.pred))
    mls.b8  <- mean(correlationVecchia::logscore(pre.b8$means, pre.b8$vars, z.pred))
    mls.cc  <- mean(correlationVecchia::logscore(pre.cc$means, pre.cc$vars, z.pred))
    
    return(list(mspe = c(mspe.b2, mspe.b3, mspe.b5, mspe.b7, mspe.b8, mspe.cc), mls = c(mls.b2, mls.b3, mls.b5, mls.b7, mls.b8, mls.cc)))
  }
}

#########################################################################

output.marginal.m.varying <- list()
for(i in 1:length(ms)) {
  output.marginal.m.varying[[i]] <- comparison(ms[i], joint = FALSE, idx = NULL)
  print(i)
}

output.joint.m.varying <- list()
for(i in 1:length(ms)) {
  output.joint.m.varying[[i]] <- comparison(ms[i], joint = TRUE, idx = NULL)
  print(i)
}

output.marginal.m.fixed <- list()
for(i in 1:length(ms)) {
  output.marginal.m.fixed[[i]] <- comparison(ms[i], joint = FALSE, idx = 10)
  print(i)
}

output.joint.m.fixed <- list()
for(i in 1:length(ms)) {
  output.joint.m.fixed[[i]] <- comparison(ms[i], joint = TRUE, idx = 10)
  print(i)
}

#########################################################################

output.marginal <- output.marginal.m.varying ; output.joint <- output.joint.m.varying
mspe.b2.marg <- c() ; mspe.b3.marg <- c() ; mspe.b5.marg <- c() ; mspe.b7.marg <- c() ; mspe.b8.marg <- c() ; mspe.cc.marg <- c()
mspe.b2.join <- c() ; mspe.b3.join <- c() ; mspe.b5.join <- c() ; mspe.b7.join <- c() ; mspe.b8.join <- c() ; mspe.cc.join <- c()
mls.b2.marg <- c() ; mls.b3.marg <- c() ; mls.b5.marg <- c() ; mls.b7.marg <- c() ; mls.b8.marg <- c() ; mls.cc.marg <- c()
for(i in 1:length(ms)) {
  
  mspe.b2.marg[i]  <- output.marginal[[i]]$mspe[1]
  mspe.b3.marg[i]  <- output.marginal[[i]]$mspe[2]
  mspe.b5.marg[i]  <- output.marginal[[i]]$mspe[3]
  mspe.b7.marg[i]  <- output.marginal[[i]]$mspe[4]
  mspe.b8.marg[i]  <- output.marginal[[i]]$mspe[5]
  mspe.cc.marg[i]  <- output.marginal[[i]]$mspe[6]
  
  mls.b2.marg[i]   <- output.marginal[[i]]$mls[1]
  mls.b3.marg[i]   <- output.marginal[[i]]$mls[2]
  mls.b5.marg[i]   <- output.marginal[[i]]$mls[3]
  mls.b7.marg[i]   <- output.marginal[[i]]$mls[4]
  mls.b8.marg[i]   <- output.marginal[[i]]$mls[5]
  mls.cc.marg[i]   <- output.marginal[[i]]$mls[6]
  
  mspe.b2.join[i]  <- output.joint[[i]]$mspe[1]
  mspe.b3.join[i]  <- output.joint[[i]]$mspe[2]
  mspe.b5.join[i]  <- output.joint[[i]]$mspe[3]
  mspe.b7.join[i]  <- output.joint[[i]]$mspe[4]
  mspe.b8.join[i]  <- output.joint[[i]]$mspe[5]
  mspe.cc.join[i]  <- output.joint[[i]]$mspe[6]
}

mspe.m.varying <- list(marginal = list(b2 = mspe.b2.marg, b3 = mspe.b3.marg, b5 = mspe.b5.marg, b7 = mspe.b7.marg, b8 = mspe.b8.marg, cc = mspe.cc.marg), 
                       joint = list(b2 = mspe.b2.join, b3 = mspe.b3.join, b5 = mspe.b5.join, b7 = mspe.b7.join, b8 = mspe.b8.join, cc = mspe.cc.join))
mls.m.varying <- list(marginal = list(b2 = mls.b2.marg, b3 = mls.b3.marg, b5 = mls.b5.marg, b7 = mls.b7.marg, b8 = mls.b8.marg, cc = mls.cc.marg), joint = NA)

#####

par(mfrow = c(1, 3))

plot(ms, sqrt(mspe.b2.join), ylim = sqrt(range(c(mspe.b5.join, mspe.cc.join))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "RMSPE", main = "RMSPE (joint)")
points(ms, sqrt(mspe.b3.join), type = "o", col = c("#984EA3"), pch = 8)
points(ms, sqrt(mspe.b5.join), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, sqrt(mspe.b7.join), type = "o", col = c("#377EB8"), pch = 15)
points(ms, sqrt(mspe.b8.join), type = "o", col = c("#FFFF33"), pch = 13)
points(ms, sqrt(mspe.cc.join), type = "o", col = c("#E41A1C"), pch = 18)

plot(ms, sqrt(mspe.b2.marg), ylim = sqrt(range(c(mspe.b5.marg, mspe.cc.marg))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "RMSPE", main = "RMSPE (marginal)")
points(ms, sqrt(mspe.b3.marg), type = "o", col = c("#984EA3"), pch = 8)
points(ms, sqrt(mspe.b5.marg), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, sqrt(mspe.b7.marg), type = "o", col = c("#377EB8"), pch = 15)
points(ms, sqrt(mspe.b8.marg), type = "o", col = c("#FFFF33"), pch = 13)
points(ms, sqrt(mspe.cc.marg), type = "o", col = c("#E41A1C"), pch = 18)

ftn <- function(x) log10(-x)

plot(ms, ftn(mls.b2.marg), ylim = ftn(range(c(mls.b5.marg, mls.cc.marg))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "log10(logscore)", main = "Log-scale logscore (marginal)")
points(ms, ftn(mls.b3.marg), type = "o", col = c("#984EA3"), pch = 8)
points(ms, ftn(mls.b5.marg), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, ftn(mls.b7.marg), type = "o", col = c("#377EB8"), pch = 15)
points(ms, ftn(mls.b8.marg), type = "o", col = c("#FFFF33"), pch = 13)
points(ms, ftn(mls.cc.marg), type = "o", col = c("#E41A1C"), pch = 18)

par(mfrow = c(1, 1))

#########################################################################

output.marginal <- output.marginal.m.fixed ; output.joint <- output.joint.m.fixed
mspe.b2.marg <- c() ; mspe.b3.marg <- c() ; mspe.b5.marg <- c() ; mspe.b7.marg <- c() ; mspe.b8.marg <- c() ; mspe.cc.marg <- c()
mspe.b2.join <- c() ; mspe.b3.join <- c() ; mspe.b5.join <- c() ; mspe.b7.join <- c() ; mspe.b8.join <- c() ; mspe.cc.join <- c()
mls.b2.marg <- c() ; mls.b3.marg <- c() ; mls.b5.marg <- c() ; mls.b7.marg <- c() ; mls.b8.marg <- c() ; mls.cc.marg <- c()
for(i in 1:length(ms)) {
  
  mspe.b2.marg[i]  <- output.marginal[[i]]$mspe[1]
  mspe.b3.marg[i]  <- output.marginal[[i]]$mspe[2]
  mspe.b5.marg[i]  <- output.marginal[[i]]$mspe[3]
  mspe.b7.marg[i]  <- output.marginal[[i]]$mspe[4]
  mspe.b8.marg[i]  <- output.marginal[[i]]$mspe[5]
  mspe.cc.marg[i]  <- output.marginal[[i]]$mspe[6]
  
  mls.b2.marg[i]   <- output.marginal[[i]]$mls[1]
  mls.b3.marg[i]   <- output.marginal[[i]]$mls[2]
  mls.b5.marg[i]   <- output.marginal[[i]]$mls[3]
  mls.b7.marg[i]   <- output.marginal[[i]]$mls[4]
  mls.b8.marg[i]   <- output.marginal[[i]]$mls[5]
  mls.cc.marg[i]   <- output.marginal[[i]]$mls[6]
  
  mspe.b2.join[i]  <- output.joint[[i]]$mspe[1]
  mspe.b3.join[i]  <- output.joint[[i]]$mspe[2]
  mspe.b5.join[i]  <- output.joint[[i]]$mspe[3]
  mspe.b7.join[i]  <- output.joint[[i]]$mspe[4]
  mspe.b8.join[i]  <- output.joint[[i]]$mspe[5]
  mspe.cc.join[i]  <- output.joint[[i]]$mspe[6]
}

mspe.m.fixed <- list(marginal = list(b2 = mspe.b2.marg, b3 = mspe.b3.marg, b5 = mspe.b5.marg, b7 = mspe.b7.marg, b8 = mspe.b8.marg, cc = mspe.cc.marg), 
                       joint = list(b2 = mspe.b2.join, b3 = mspe.b3.join, b5 = mspe.b5.join, b7 = mspe.b7.join, b8 = mspe.b8.join, cc = mspe.cc.join))
mls.m.fixed <- list(marginal = list(b2 = mls.b2.marg, b3 = mls.b3.marg, b5 = mls.b5.marg, b7 = mls.b7.marg, b8 = mls.b8.marg, cc = mls.cc.marg), joint = NA)

#####

par(mfrow = c(1, 3))

plot(ms, sqrt(mspe.b2.join), ylim = sqrt(range(c(mspe.b5.join, mspe.cc.join))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "RMSPE", main = "RMSPE (joint)")
points(ms, sqrt(mspe.b3.join), type = "o", col = c("#984EA3"), pch = 8)
points(ms, sqrt(mspe.b5.join), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, sqrt(mspe.b7.join), type = "o", col = c("#377EB8"), pch = 15)
points(ms, sqrt(mspe.b8.join), type = "o", col = c("#FFFF33"), pch = 13)
points(ms, sqrt(mspe.cc.join), type = "o", col = c("#E41A1C"), pch = 18)

plot(ms, sqrt(mspe.b2.marg), ylim = sqrt(range(c(mspe.b5.marg, mspe.cc.marg))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "RMSPE", main = "RMSPE (marginal)")
points(ms, sqrt(mspe.b3.marg), type = "o", col = c("#984EA3"), pch = 8)
points(ms, sqrt(mspe.b5.marg), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, sqrt(mspe.b7.marg), type = "o", col = c("#377EB8"), pch = 15)
points(ms, sqrt(mspe.b8.marg), type = "o", col = c("#FFFF33"), pch = 13)
points(ms, sqrt(mspe.cc.marg), type = "o", col = c("#E41A1C"), pch = 18)

ftn <- function(x) log10(-x)

plot(ms, ftn(mls.b2.marg), ylim = ftn(range(c(mls.b5.marg, mls.cc.marg))), type = "o", col = c("#4DAF4A"), pch = 16, xlab = "increasing m", ylab = "log10(logscore)", main = "Log-scale logscore (marginal)")
points(ms, ftn(mls.b3.marg), type = "o", col = c("#984EA3"), pch = 8)
points(ms, ftn(mls.b5.marg), type = "o", col = c("#FF7F00"), pch = 17)
points(ms, ftn(mls.b7.marg), type = "o", col = c("#377EB8"), pch = 15)
points(ms, ftn(mls.b8.marg), type = "o", col = c("#FFFF33"), pch = 13)
points(ms, ftn(mls.cc.marg), type = "o", col = c("#E41A1C"), pch = 18)

par(mfrow = c(1, 1))

#########################################################################

save(locs, z, X.train, locs.pred, z.pred, X.pred, ns_obs, ns_pred,
     ms, output.joint.m.varying, output.joint.m.fixed, output.marginal.m.varying, output.marginal.m.fixed, 
     mspe.m.varying, mls.m.varying, mspe.m.fixed, mls.m.fixed, table.runtime,
     file = "DATA/joint_CRCM_NCEP_10152021_prediction.RData")
