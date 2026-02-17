library(tidyverse)

sup_norm_diff <- function(step_approx, times, true_baseline=function(x){x}){
    # for step function against a constant baseline hazard
    true_comparison_points <- true_baseline(times)
    return(max(c(abs(step_approx - true_comparison_points), 
                 abs(step_approx - c(0, true_comparison_points[1:(length(times)-1)])))))
}


step_functions_sup_diff <- function(t1, v1, t2, v2){
    t1 <- c(0, t1); v1 <- c(0, v1); t2 <- c(0, t2); v2 <- c(0, v2)
    minmaxtime <- min(t1[length(t1)], t2[length(t2)])
    timepoints <- sort(unique(c(t1, t2)))
    minmax_index <- findInterval(minmaxtime, timepoints)
    timepoints <- timepoints[1:minmax_index]
    curve1 <- v1[findInterval(timepoints, t1)]
    curve2 <- v2[findInterval(timepoints, t2)]
    return(max(abs(curve1 - curve2)))
}


cdp_beta_experiment <- function(nsamples, beta, epsilon, censor=0.3, C_z=1,
                                delta=0.001, reps=10, sensitivity=NA, stepsize=0.5,
                                sens_Cz=1){
    beta_error <- rep(0, reps)
    for (r in 1:reps){
        obs_beta <- sim_observations(nsamples, beta, censor=censor, C_z=C_z)
        private_beta <- cdp_cox(obs_beta, epsilon, delta=delta, C_z=sens_Cz, 
                                sensitivity=sensitivity, stepsize=stepsize)
        beta_error[r] <- sum((beta-private_beta)^2)
    }
    return(beta_error)
}


cdp_full_experiment <- function(nsamples, beta, epsilon, censor=0.3, C_z=1,
                                delta=0.001, reps=10, true_Lambda=function(x){x},
                                max_time=1, split=c(1, 0.1, 1), sensitivity=NA,
                                truncation=NA){
    beta_error <- rep(0, reps)
    lambda_error <- rep(0, reps)
    for (r in 1:reps){
        obs_beta <- sim_observations(as.integer(nsamples*split[1]), beta, 
                                     censor=censor, C_z=C_z)
        private_beta <- cdp_cox(obs_beta, epsilon, delta=delta, C_z=C_z, 
                                sensitivity=sensitivity)
        
        obs_p <- sim_observations(as.integer(nsamples*split[2]), beta, 
                                  censor=censor, C_z=C_z)
        
        p_hat <- fdp_probabilities(list(obs_p[, 1]), c(epsilon))

        obs_tree <- sim_observations(as.integer(nsamples*split[3]), 
                                     beta, censor=censor, C_z=C_z)
        
        breslow <- fdp_breslow(list(obs_tree), c(epsilon), private_beta, p_hat, 
                               C_z=C_z, delta=c(delta), cutoff=max_time,
                               truncation=truncation)
        beta_error[r] <- sum((private_beta - beta)^2)
        lambda_error[r] <- sup_norm_diff(breslow$vals, breslow$times,
                                         true_baseline=true_Lambda) 
    }
    return(list('beta_error'=beta_error, 'lambda_error'=lambda_error))
}


fdp_experiment <- function(S, nsamples, beta, epsilon, censor=0.3, C_z=1, C_beta=1,
                           stepsize=0.5, delta=rep(0.001, S), reps=10, max_time=1, 
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
            priv_beta <- fdp_cox_interactive(obs, epsilon, delta=delta, 
                                             stepsize=stepsize, C_z=C_z)
        } else{
            priv_beta <- fdp_cox(obs, epsilon, stepsize=stepsize, delta=delta, 
                                 C_z=C_z)
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
                               cutoff=max_time, C_z=C_z)
        sup_errors[r] <- sup_norm_diff(priv_breslow$vals, priv_breslow$times,
                                       true_baseline=true_Lambda)
    }
    return(list('beta_error'=beta_errors, 'lambda_error'=sup_errors))
}



line_plot <- function(results, x_axis, factors, xlab='Samples', 
                      colour=expression(epsilon), ylab='Squared error', 
                      bar=1500, factornames=NA, ymax=NA, xmax=NA, df_func=F){
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

    if(df_func){
        return(df)
    }
    
    ggplot(df, aes(x=x, y=mean, color=s, group=s)) +
        geom_line(size = 0.9) +
        geom_point(size = 2) +
        scale_x_continuous(breaks=x_axis) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=bar, linewidth=0.7,
                      alpha=0.6) +
        labs(x=xlab, y=ylab, colour=colour) + 
        theme_minimal() +
        coord_cartesian(ylim=c(min(0, df$mean), max(ymax, max(df$upper))))
}
