source('experiment_functions.R')

set.seed(123)
reps <- 200

## CDP experiments 
epsilons <- c(0.75, 1:6)
nsamples <- seq(from=20000, to=50000, by=5000)

cdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
cdp_breslow <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps)
        cdp_beta[e, s, ] <- results$beta_error
        cdp_breslow[e, s, ] <- results$lambda_error
    }
}

# varying dimension
param_dim <- c(2, 3, 4, 5, 6, 7, 8)
cdp_dimension <- array(0, dim=c(length(param_dim), length(epsilons), reps))
for (e in 1:length(param_dim)){
    print(param_dim[e])
    for (s in 1:length(epsilons)){
        print(epsilons[s])
        for (r in 1:reps){
            beta <- rnorm(param_dim[e])
            if (sum(beta^2) > 1){
                beta <- 1/sqrt(sum(beta^2)) * beta
            }
            obs <- sim_observations(30000, beta, censor=0.3)
            out <- cdp_sgd(obs, epsilons[s])
            cdp_dimension[e, s, r] <- sum((beta - out)^2)
        }
    }
}
cdp_dimension <- aperm(cdp_dimension, perm = c(2, 1, 3))


# sensitivity analysis
sens <- seq(from=0.005, to=0.02, by=0.0025)
sens_sensitivity <- array(0, dim=c(length(epsilons), length(sens), reps))
for (e in 1:length(epsilons)){
    for (s in 1:length(sens)){
        sens_sensitivity[e, s, ] <- cdp_beta_experiment(30000, c(0, 0.5, 0.8),
                                                 epsilon=epsilons[e], reps=reps,
                                                 sensitivity=sens[s])
    }
}

step <- seq(0.2, 0.8, by=0.1)
step_sensitivity <- array(0, dim=c(length(epsilons), length(step), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(step)){
        print(cdp_beta_experiment(30000, c(0, 0.5, 0.8),
                                  epsilon=epsilons[e], reps=reps,
                                  stepsize=step[s]))
        step_sensitivity[e, s, ] <- cdp_beta_experiment(30000, c(0, 0.5, 0.8),
                                                        epsilon=epsilons[e], reps=reps,
                                                        stepsize=step[s])
    }
}

# varying censoring rate 
censoring_rate <- seq(from=0.1, to=1.3, by=0.2)
censoring_beta <- array(0, dim=c(length(epsilons), length(censoring_rate), reps))
censoring_breslow <- array(0, dim=c(length(epsilons), length(censoring_rate), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(censoring_rate)){
        print(censoring_rate[s])
        results <- cdp_full_experiment(30000, c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, censor=censoring_rate[s])
        censoring_beta[e, s, ] <- results$beta_error
        censoring_breslow[e, s, ] <- results$lambda_error
    }
}

## FDP experiments
nservers1 <- 2:8
# experiments where data is reused for fully-interactive mechanism. 
fdp_beta <- array(0, dim=c(length(epsilons), length(nservers), reps))
fdp_hazard <- array(0, dim=c(length(epsilons), length(nservers), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nservers)){
        print(nservers[s])
        results <- fdp_experiment(nservers1[s], rep(10000, nservers[s]), c(0, 0.5, 0.8), 
                                  rep(epsilons[e], nservers[s]), censor=0.3, reps=reps,
                                  interative=TRUE)
        fdp_beta[e, s, ] <- results$beta_error
        fdp_hazard[e, s, ] <- results$lambda_error
    }
}

nservers2 <- c(2, 4, 8, 12, 16, 20)
# experiments satisfying the FDP constraints 
fdps_beta <- array(0, dim=c(length(epsilons), length(nservers), reps))
fdps_hazard <- array(0, dim=c(length(epsilons), length(nservers), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nservers)){
        print(nservers[s])
        results <- fdp_experiment(nservers2[s], rep(10000, nservers[s]), c(0, 0.5, 0.8), 
                                  rep(epsilons[e], nservers[s]), censor=0.3, 
                                  reps=reps, interactive=FALSE)
        fdps_beta[e, s, ] <- results$beta_error
        fdps_hazard[e, s, ] <- results$lambda_error
    }
}

# code for plotting
line_plot(cdp_beta, nsamples, epsilons, title='Coefficient estimation: varying no. of samples')
ggsave('cdp_beta.png', width=5, height=4)
line_plot(cdp_breslow, nsamples, epsilons, ylab='Sup error', title='Cumulative hazard: varying no. of samples')
ggsave('cdp_breslow.png', width=5, height=4)
line_plot(cdp_dimension, param_dim, epsilons, bar=0.2, 
          xlab='Dimension', title='Coefficient estimation: varying dimension')
ggsave('cdp_dimension.png', width=6, height=4)

line_plot(sens_sensitivity, sens, epsilons, bar=0.001, 
          xlab='Sensitivity scaling')
ggsave('noise_sensitivity.png', width=6, height=4)

line_plot(step_sensitivity, step, epsilons, bar=0.05, xlab='Step size')
ggsave('stepsize_sensitivity.png', width=6, height=4)

line_plot(censoring_beta, censoring_rate, epsilons, bar=0.05,
          xlab='Censoring distribution rate')
ggsave('censoring_beta.png', width=6, height=4)
line_plot(censoring_breslow, censoring_rate, epsilons, bar=0.05,
          xlab='Censoring distribution rate', ylab='Sup error')
ggsave('censoring_breslow.png', width=6, height=4)


line_plot(fdp_beta, nservers, epsilons, bar=0.3, xlab='Servers')
ggsave('fdpreuse_beta.png', width=6, height=4)
line_plot(fdp_hazard, nservers, epsilons, ylab='Sup error', xlab='Servers', bar=0.3)
ggsave('fdpreuse_breslow.png', width=6, height=4)
line_plot(fdps_beta, nservers, epsilons, bar=0.3, xlab='Servers')
ggsave('fdpalg_beta.png', width=6, height=4)
line_plot(fdps_hazard, nservers, epsilons, ylab='Sup error', xlab='Servers', bar=0.3)
ggsave('fdpalg_breslow.png', width=6, height=4)

line_plot(fdps_beta, nservers, epsilons, bar=0.3, xlab='Servers')
ggsave('fdpalg_beta.png', width=6, height=4)
line_plot(fdps_hazard, nservers, epsilons, ylab='Sup error', xlab='Servers', bar=0.3)
ggsave('fdpalg_breslow.png', width=6, height=4)

line_plot(fdplarge_beta, c(2, 4*(1:5)), epsilons, bar=1, xlab='Servers')
ggsave('fdp25_beta.png', width=6, height=4)
line_plot(fdplarge_hazard, c(2, 4*(1:5)), epsilons, bar=1, xlab='Servers',
          ylab='Sup error')
ggsave('fdp25_hazard.png', width=6, height=4)



#testing refactor 
result <- fdp_experiment(4, rep(100000, 4), c(0, 0.5, 0.8), rep(100, 4),
                         censor=0.3, interactive=FALSE)
