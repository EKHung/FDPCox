library(tidyverse)

fdp_cox_interactive <- function(obs, epsilon, delta=0.001, niters=NA, C_beta = 1,
                                C_z = 1, sensitivity=NA, stepsize=0.5, 
                                weights=NA){
    # input obs is a list indexed [[s]][i, (t, \delta, preds)]
    # requires number of preds to be the same across the different servers
    # epsilon, sensitivity, and weights are length S positive vectors
    
    S <- length(obs); p <- dim(obs[[1]])[2] - 2;
    nsamples <- rep(0, S)
    sorted_obs <- list()
    for (s in 1:S){
        nsamples[s] <- dim(obs[[s]])[1]
        sorted_obs[[s]]  <- preprocess(obs[[s]])
    }
    
    if (is.na(niters)){
        niters <- 6*log(C_z^2 * sum(nsamples)/p^2)
    }
    if (is.na(weights)){
        weights <- pmin(nsamples, (nsamples^2*epsilon^2/p))
        weights <- weights/sum(weights)
    }
    if (is.na(sensitivity)){
        sensitivity <- sapply(nsamples, function(n){
            4*C_z/n + 5*exp(2*C_z*C_beta)*max(C_z, C_z^2)*log(n+1)/n})
    }
    
    noise_scale <- sqrt(sensitivity^2*niters/epsilon * (2*log(1/delta)/epsilon + 1))
    
    beta <- rep(0, p)
    for (k in 1:niters){
        step <- 0
        for (s in 1:S){
            step <- step + weights[s]*((gradient(beta, sorted_obs[[s]]))
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

sup_norm_diff <- function(step_approx, times, true_baseline=function(x){x}){
    # for step function against a constant baseline hazard
    true_comparison_points <- true_baseline(times)
    return(max(c(abs(step_approx - true_comparison_points), 
                 abs(step_approx - c(0, true_comparison_points[1:(length(times)-1)])))))
}


cdp_beta_experiment <- function(nsamples, beta, epsilon, censor=0.3, C_z=1,
                                delta=0.001, reps=10, sensitivity=NA, stepsize=0.5){
    beta_error <- rep(0, reps)
    for (r in 1:reps){
        #cat(paste("\rrepetition", r))
        obs_beta <- sim_observations(nsamples, beta, censor=censor, C_z=C_z)
        private_beta <- cdp_cox(obs_beta, epsilon, delta=delta, C_z=C_z, 
                                sensitivity=sensitivity, stepsize=stepsize)
        beta_error[r] <- sum((beta-private_beta)^2)
    }
    return(beta_error)
}


cdp_full_experiment <- function(nsamples, beta, epsilon, censor=0.3, C_z=1,
                                delta=0.001, reps=10, true_Lambda=function(x){x},
                                max_time=1, split=c(1, 0.1, 1)){
    beta_error <- rep(0, reps)
    lambda_error <- rep(0, reps)
    for (r in 1:reps){
        cat(paste("\rrepetition", r))
        obs_beta <- sim_observations(nsamples*split[1], beta, censor=censor, C_z=C_z)
        private_beta <- cdp_cox(obs_beta, epsilon, delta=delta, C_z=C_z)
        
        obs_p <- sim_observations(nsamples*split[2], beta, censor=censor, C_z=C_z)
        p_hat <- fdp_probabilities(list(obs_p), c(epsilon))

        obs_tree <- sim_observations(nsamples*split[3], beta, censor=censor, C_z=C_z)
        breslow <- fdp_breslow(list(obs_tree), c(epsilon), private_beta, p_hat, 
                               C_z=C_z, delta=c(delta), cutoff=max_time)
        beta_error[r] <- sum((private_beta - beta)^2)
        lambda_error[r] <- sup_norm_diff(breslow$vals, breslow$times,
                                         true_baseline=true_Lambda) 
    }
    return(list('beta_error'=beta_error, 'lambda_error'=lambda_error))
}


fdp_experiment <- function(S, nsamples, beta, epsilon, censor=1, C_z=1, C_beta=1,
                           stepsize=0.5, delta=NA, reps=10, max_time=1, 
                           split=c(1, 0.1, 1), true_Lambda=function(x){x},
                           interactive=FALSE){
    # setting interactive=TRUE uses the whole set of data at each gradient descent
    # iteration to make use of RDP composition. 
    
    beta_errors <- rep(0, reps)
    sup_errors <- rep(0, reps)
    
    for (r in 1:reps){
        #cat(paste("\rrepetition", r))
        obs <- list()
        # beta estimation part 
        for (s in 1:S){
            obs[[s]] <- sim_observations(nsamples[s]*split[1], beta, censor=censor, 
                                         C_z=C_z)
        }
        if (interactive){
            priv_beta <- fdp_cox_interactive(obs, epsilon, stepsize=stepsize, C_z=C_z)
        } else{
            priv_beta <- fdp_cox(obs, epsilon, stepsize=stepsize, C_z=C_z)
        }
        beta_errors[r] <- sum((beta - priv_beta)^2)
        
        prob_obs <- list()
        # at-risk prob part 
        for (s in 1:S){
            prob_obs[[s]] <- sim_observations(nsamples[s]*split[2], beta, 
                                              censor=censor, C_z=C_z)[, 1]
        }
        p_hat <- fdp_probabilities(obs, epsilon, delta=delta, cutoff=max_time)
        
        # cumulative_hazard estimation part
        lambda_obs <- list()
        for (s in 1:S){
            lambda_obs[[s]] <- sim_observations(nsamples[s]*split[3], beta, 
                                                censor=censor, C_z=C_z)
        }
        priv_breslow <- fdp_breslow(obs, epsilon, priv_beta, p_hat, delta=delta, 
                               cutofff=max_time, C_z=C_z)
        sup_errors[r] <- sup_norm_diff(priv_breslow$vals, priv_breslow$times,
                                       true_baseline=true_Lambda)
    }
    return(list('beta_error'=beta_errors, 'lambda_error'=sup_errors))
}



line_plot <- function(results, x_axis, factors, xlab='Samples', 
                      colour=expression(epsilon), ylab='Squared error', bar=1500){
    results_mean <- t(apply(results, c(1, 2), mean))
    results_sd <- t(apply(results, c(1, 2), sd))
    results_sd <- results_sd / sqrt(dim(results)[3])
    results_lower <- results_mean - results_sd
    results_upper <- results_mean + results_sd
    df <- data.frame(
        x = rep(x_axis, times = length(factors)),
        s = factor(rep(factors, each = length(x_axis))),
        mean = as.vector(results_mean),
        lower = as.vector(results_lower),
        upper = as.vector(results_upper)
    )
    
    ggplot(df, aes(x=x, y=mean, color=s, group=s)) +
        geom_line(size = 0.9) +
        geom_point(size = 1.8) +
        scale_x_continuous(breaks=x_axis) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=bar, alpha=0.4) +
        ylim(0, max(df[, c('mean', 'lower', 'upper')])) +
        if (colour == 'eps'){
            labs(x=xlab, y=ylab, color=expression(epsilon))
        } else{
            labs(x=xlab, y=ylab, color=colour)
        }
}
