sim_CoxPH <- function(beta, covariates){
    # Assumes that the baseline hazard rate $\lambda_0(t) = 1$
    u <- runif(1)
    return(-log(u)*exp(-sum(beta*covariates)))
}


preprocess <- function(obs, cutoff=1){
    # obs is a n * (2 + p) matrix with columns (T, \Delta, Z)
    # Return a list containing:
        # obs = n * (2+p) matrix sorted by observation times
        # uncensored = indices of obs with \delta_i =1
        # u_preds = 1/n \sum_{\delta_i = 1} Z_i
    p <- dim(obs)[2]
    sorted_obs <- obs[order(obs[, 1]), ]
    if (is.na(cutoff)){
      uncensored <- which(sorted_obs[, 2] == 1)
    } else{
      cutoff_ind <- findInterval(cutoff, sorted_obs[, 1])
      uncensored <- which(sorted_obs[1:cutoff_ind, 2] == 1)
    }
    uncensored_preds <- apply(sorted_obs[uncensored, 3:p], MARGIN=2, sum)
    return(list(obs=sorted_obs, uncensored=uncensored, upreds=uncensored_preds))
}


gradient <- function(beta, obs){
    # obs is a list output from preprocess()
    # Return gradient of normalised log-likelihood from obs evaluated at beta
    sorted_obs <- obs$obs
    n <- nrow(sorted_obs); p <- ncol(sorted_obs)
    weights <- exp(sorted_obs[, 3:p] %*% beta)
    tsum_weights <- rev(cumsum(rev(weights)))
    weighted_preds <- sweep(sorted_obs[, 3:p], MARGIN=1, weights, '*')
    tsum_preds <- apply(weighted_preds, 2, function(x) rev(cumsum(rev(x))))
    grad <- obs$upreds
    for (t in obs$uncensored){
        grad <- grad - tsum_preds[t, ] / tsum_weights[t]
    }
    return(grad/n)
}


cdp_tree <- function(data, beta_hat, epsilon, levels, truncation, delta=0.001, 
                     max_time=1){
  sorted_obs <- preprocess(data)
  n <- dim(data)[1]; p <- dim(data)[2]
  weights <- exp(sorted_obs$obs[, 3:p] %*% beta_hat)
  tsum_weights <- rev(cumsum(rev(weights)))
  tsum_weights <- pmin(1/tsum_weights, 1/(truncation*n))

  # non-discretised breslow estimator at beta_hat
  timepoints <- length(sorted_obs$uncensored)
  breslow_beta <- list('time'=rep(0, timepoints), 'breslow'=rep(0, timepoints))
  for (u in 1:timepoints){
    breslow_beta$time[u] <- sorted_obs$obs[sorted_obs$uncensored[u], 1]
    breslow_beta$breslow[u] <- tsum_weights[sorted_obs$uncensored[u]]
  }
  
  tree <- matrix(0, nrow=levels, ncol=2^levels)
  counter <- 1; old_counter <- 1
  for (m in 1:(2^levels)){
    while(breslow_beta$time[counter]<m/(2^levels) && counter<timepoints){
      counter <- counter + 1
    }
    tree[levels, m] <- sum(breslow_beta$breslow[old_counter:(counter-1)])
    old_counter <- counter
  }
  
  for (l in (levels-1):1){
    tree_level <- rep(0, 2^(levels-l))
    for (j in 1:2^l){
      tree[l, j] <- tree[l+1, 2*j-1] + tree[l+1, 2*j]
    }
  }
  noise_scale <- ((2*log(1/delta)/epsilon + 1)/(epsilon*n^2/levels)
                  *(1/truncation^4 + 3/truncation^2))
  for(l in 1:levels){
    for(m in 1:2^l){
      tree[l, m] <- tree[l, m] + rnorm(1, sd=sqrt(noise_scale))
    }
  }
  return(tree)
}


fdp_trees_to_breslow <- function(trees_list, weights){
  S <- length(trees_list)
  height <- dim(trees_list[[1]])[1]; bins <- dim(trees_list[[1]])[2]
  timepoints <- (1:bins)/bins
  steps <- rep(0, bins)
  for(t in 1:bins){
    binary_rep <- rev(as.integer(intToBits(t))[1:height]) 
    for (l in 1:height){
      if (binary_rep[l] == 1){
        for (s in 1:S){
          steps[t] <- steps[t] + 
            weights[s]*trees_list[[s]][l, sum(2^((l:1) - 1)*binary_rep[1:l])]
        }
      }
    }
  }
  for (s in 1:S){
    steps[bins] <- steps[bins] + weights[s]*sum(trees_list[[s]][1, ])
  }
  return(list('times'=timepoints, 'vals'=steps))
}



