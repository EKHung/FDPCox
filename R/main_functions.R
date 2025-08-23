#' Generate a dataset
#'
#' Generates observations, with i.i.d. covariates from a uniform distribution, 
#' and a failure time from the Cox proportional hazards model with baseline hazard 1.
#' The recorded observation may be censored according to an exponential distribution. 
#'
#' @param n Integer. Number of observations to simulate.
#' @param beta Numeric vector. d-dimensional vector of the Cox coefficients.
#' @param censor Optional positive numeric. Censoring distribution is Exp(censor). Defaults to 1.
#' @param C_z Optional positive numeric. Covariates for each observation are generated as 
#' i.i.d. from Uniform(-C_z/sqrt(d), C_z/sqrt(d)). Defaults to 1.
#' @return A n*(2+length(beta)) matrix, where the first column is the observed time,
#' the second column is the 1/0 (failure observed/censored) and remaining columns are covariates.
#' @export
sim_observations <- function(n, beta, censor=1, C_z=1){
    # censoring distribution is Exp(censor)
    d <- length(beta)
    predictors <- matrix(runif(n*d, min=-C_z/sqrt(d), max=C_z/sqrt(d)),
                         nrow=n, ncol=d)
    
    obs <- matrix(0, nrow=n, ncol=4)
    obs[, 1] <- rexp(n, censor)
    obs[, 2] <- apply(predictors, 1, 
                      function(x) return(rexp(1, exp(sum(beta*x)))))
    obs[, 3] <- apply(obs[, c(1, 2)], 1, min)
    obs[, 4] <- obs[, 1] > obs[, 2] # indicators for not censored
    return(cbind(obs[, c(3, 4)], predictors))
}


#' CDP Cox regression
#'
#' A private SGD implementation to estimate Cox regression coefficients, that satisfies
#' central differential privacy. 
#'
#' @param obs n*(2+d) dimensional matrix. Output from sim_observations(). 
#' @param epsilon Positive numeric. Epsilon privacy budget in terms of DP.
#' @param delta Optional positive numeric. Delta privacy budget in terms of DP. 
#' Defaults to 0.001. 
#' @param niters Optional positive integer. Number of gradient descent iterations. 
#' Defaults to 6*log(C_z^2 *n/p^2). 
#' @param C_beta Optional positive numeric. l_2 upper bound for coefficients. Defaults to 1. 
#' @param C_z Optional positive numeric. l_2 upper bound for covariates. Defaults to 1. 
#' @param sensitivity Optional positive numeric. l_2-sensitivity of the gradient of the partial 
#' log-likelihood. Defaults to an upper bound for ensuring privacy based on the other parameter values. 
#' @param stepsize Optional positive numeric. Step size in the descent direction. Defaults to 0.5. 
#' @return A d-dimensional numeric vector estimating the Cox regression coefficients. 
#' @export
cdp_cox <- function(obs, epsilon, delta=0.001, niters=NA, C_beta = 1, C_z = 1, 
                    sensitivity=NA, stepsize=0.5, cutoff=1){
    sorted_obs <- preprocess(obs, cutoff=cutoff)
    n <- dim(obs)[1]; p <- dim(obs)[2] - 2
    
    if (is.na(niters)){
        niters <- 6*log(C_z^2 *n/p^2)
    }
    if (is.na(sensitivity)){
        sensitivity <- 4*C_z/n + 5*exp(2*C_z*C_beta)*max(C_z, C_z^2)*log(n+1)/n
    }
    noise_scale <- sqrt(sensitivity^2*niters/epsilon * (2*log(1/delta)/epsilon + 1))
    beta <- rep(0, p)
    for (k in 1:niters){
        beta <- beta + stepsize*(gradient(beta, sorted_obs) + noise_scale * rnorm(p))
        # projection onto l_2 ball of radius C_beta
        if (sqrt(sum(beta^2)) > C_beta){
            beta <- C_beta / sqrt(sum(beta^2)) * beta
        }
    }
    return(beta)
}


#' FDP marginal survival estimation
#'
#' Marginal survival estimation via mean estimation with federated differential 
#' privacy guarantees. 
#'
#' @param times list of S non-negative numeric vectors.  
#' @param epsilon Positive numeric vector of length S. Epsilon privacy budget in terms of DP.
#' @param delta Optional positive numeric vector of length S. Delta privacy budget 
#' in terms of DP. Defaults to S entries of 0.001.
#' @param cutoff. Optional positive numeric. Timepoint at which to estimate probability 
#' of being at-risk. Defaults to 1. 
#' @return A numeric that is an FDP estimate of P(Y(1)=cutoff). 
#' @export
fdp_probabilities <- function(times, epsilon, delta=NA, cutoff=1){
    S <- length(times)
    if (is.na(delta)){
        delta <- rep(0.001, S)
    }
    p0 <- 0
    for (s in 1:S){
        noise <- 2*log(1.25/delta[s])/(length(times[[s]])^2 * epsilon[s]^2)
        p0 <- p0 + mean(times[[s]]>1) + rnorm(1, sd=sqrt(noise))
    }
    return(p0)
}


#' FDP Breslow
#'
#' A distributed variant of the Breslow estimator satisfying federated differential privacy. 
#'
#' @param datasets List of S datasets, with (2+d) columns, as generated from sim_observations(). 
#' @param epsilon Positive numeric vector of length S. Epsilon privacy budget in terms of DP. 
#' @param beta_hat Numeric vector of length d. Estimate of the Cox regression coefficients 
#' @param p_hat Positive numeric. Estimate of P(Y(cutoff)=1). 
#' @param delta Optional positive numeric vector of length S. Delta privacy budget 
#' in terms of DP. Defaults to S entries of 0.001.
#' @param cutoff Optional positive numeric. Defaults to 1. 
#' @param weights Optional positive numeric vector of length S. Weights used to aggregate 
#' the per-server cumulative hazard functions. Defaults to effective sample sizes, normalised
#' to sum to 1. 
#' @param tree_height Optional positive numeric. Number of tree levels to 
#' @return list containing two positive numeric vectors, 'times' and 'vals'. 
#' vals[t] is the estimate for the cumulative hazard at times[t]. 
#' @export
fdp_breslow <- function(datasets, epsilon, beta_hat, p_hat, delta=NA, cutofff=1, 
                        C_z=1, weights=NA, tree_height=NA){
    S <- length(datasets)
    if (is.na(delta)){
        delta <- rep(0.001, S)
    }
    nsamples <- rep(0, S)
    for (s in 1:S){
        nsamples[s] <- dim(datasets[[s]])[1]
    }
    if (is.na(weights)){
        weights <- pmin(nsamples, nsamples^2*epsilon^2)
        weights <- weights/sum(weights)
    } 
    if (is.na(tree_height)){
        tree_height <-  floor(log(sum(pmin(nsamples, nsamples^2*epsilon^2)), 2)/2)
    }
    trees_list <- list()
    truncation <- exp(-C_z*sqrt(sum(beta_hat)^2))*p_hat*0.9

    for (s in 1:S){
        trees_list[[s]] <- cdp_tree(datasets[[s]], beta_hat, epsilon[s], truncation=truncation,
                                    levels=tree_height, delta=delta[s], max_time=cutoff)
    }
    priv_breslow <- fdp_trees_to_breslow(trees_list, weights)
    return(priv_breslow)
}


#' FDP Cox regression
#'
#' A private SGD implementation to estimate Cox regression coefficients, that satisfies
#' federated differential privacy. 
#' 
#' @param datasets List of S datasets, with (2+d) columns, as generated from sim_observations(). 
#' @param epsilon Positive numeric vector of length S. Epsilon privacy budget in terms of DP. 
#' @param delta Optional positive numeric vector of length S. Delta privacy budget 
#' in terms of DP. Defaults to S entries of 0.001.
#' @param weights. Optional positive numeric vector of length S. Weights used to aggregate 
#' the gradient estimates. Defaults to effective sample sizes, normalised to sum to 1. 
#' @param niters Optional positive integer. Number of gradient descent iterations. 
#' Defaults to 6*log(C_z^2 *sum(samples)/p^2). 
#' @param C_beta Optional positive numeric. l_2 upper bound for coefficients. Defaults to 1. 
#' @param C_z Optional positive numeric. l_2 upper bound for covariates. Defaults to 1. 
#' @param sensitivity Optional positive numeric vector of length S. l_2-sensitivity of 
#' the gradient of the partial log-likelihood. Defaults to an upper bound for ensuring
#' privacy based on the other parameter values. 
#' @param stepsize Optional positive numeric. Step size in the descent direction. Defaults to 0.5. 
#' @return A d-dimensional numeric vector estimating the Cox regression coefficients. 
#' @export
fdp_cox <- function(obs, epsilon, delta=NA, niters=NA, C_beta=1, C_z=1,
                    sensitivity=NA, stepsize=0.5, weights=NA){
  
  S <- length(obs); p <- dim(obs[[1]])[2] - 2;
  nsamples <- rep(0, S)
  for (s in 1:S){
    nsamples[s] <- dim(obs[[s]])[1]
  }
  if (is.na(delta)){
      delta <- rep(0.001, S)
  }
  if (is.na(niters)){
    niters <- 6*log(C_z^2 * sum(nsamples) / p^2)
  }
  if (is.na(weights)){
    weights <- pmin(nsamples, (nsamples^2*epsilon^2/p))
    weights <- weights/sum(weights)
  }
  nsamples <- nsamples %/% niters  #set batch size
  if (is.na(sensitivity)){
    sensitivity <- sapply(nsamples, function(n){
      4*C_z/n + 5*exp(2*C_z*C_beta)*max(C_z, C_z^2)*log(n+1)/n})
  }
  noise_scale <- sqrt(2*log(1.25/delta))*sensitivity/epsilon
  
  beta <- rep(0, p)
  for (k in 1:niters){
    step <- 0
    for (s in 1:S){
      subset_obs <- obs[[s]][(1+(k-1)*nsamples[s]):(k*nsamples[s]), ]
      sorted_obs <- preprocess(subset_obs)
      step <- step + weights[s]*((gradient(beta, sorted_obs))
                                 + noise_scale[s]*rnorm(p))
    }
    beta <- beta + stepsize*step
    # projection onto l_2 ball of radius C_beta
    if (sqrt(sum(beta^2)) > C_beta){
      beta <- C_beta / sqrt(sum(beta^2)) * beta
    }
  }
  return(beta)
}


