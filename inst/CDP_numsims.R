library('FDPCox')
source('experiment_functions.R')

set.seed(123)

reps <- 100
epsilons <- c(0.2, 1, 5, 15)
nsamples <- 2^(12:18)
lognumcdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
lognumcdp_breslow <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, sensitivity=NA)
        lognumcdp_beta[e, s, ] <- results$beta_error
        lognumcdp_breslow[e, s, ] <- results$lambda_error
    }
}

#without covariates 
lognumcdp_NA <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        for (r in 1:reps){
            obs_p <- sim_observations(as.integer(nsamples[s]*0.1), c(0), 
                                      censor=0.3, C_z=1)
            p_hat <- fdp_probabilities(list(obs_p[, 1]), c(epsilons[e]), 
                                       delta=c(0.001))
            
            obs_tree <- sim_observations(as.integer(nsamples[s]), c(0), 
                                         censor=0.3, C_z=1)
            
            breslow <- fdp_breslow(list(obs_tree), c(epsilons[e]), c(0), p_hat, 
                                   C_z=1, delta=c(0.001), cutoff=1)
            lognumcdp_NA[e, s, r] <-  sup_norm_diff(breslow$vals, breslow$times,
                                                    true_baseline=function(x){x}) 
        }
    }
}

# log plots for samples
line_plot(log(lognumcdp_beta, 2), round(log(nsamples, 2), 2), epsilons, bar=0.05,
          ylab='log(Square error)', xlab='log(samples)')
ggsave('cdp_beta_logsamples.png', width=4, height=3)
line_plot(log(lognumcdp_breslow, 2), round(log(nsamples, 2), 2), epsilons, bar=0.05,
          ylab='log(Sup error)', xlab='log(samples)')
ggsave('cdp_breslow_logsamples.png', width=4, height=3)
line_plot(log(lognumcdp_NA, 2), round(log(nsamples, 2), 2), epsilons, bar=0.05,
          ylab='log(Sup error)', xlab='log(samples)')
ggsave('cdp_NA_logsamples.png', width=4, height=3)


epsilons <- 2^seq(-4, 10, 2)
nsamples <- c(10000, 20000, 50000)
logepscdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
logepscdp_breslow <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, sensitivity=NA)
        logepscdp_beta[e, s, ] <- results$beta_error
        logepscdp_breslow[e, s, ] <- results$lambda_error
    }
}

logepscdp_NA <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        for (r in 1:reps){
            obs_p <- sim_observations(as.integer(nsamples[s]*0.1), c(0), 
                                      censor=0.3, C_z=1)
            p_hat <- fdp_probabilities(list(obs_p[, 1]), c(epsilons[e]), 
                                       delta=c(0.001))
            
            obs_tree <- sim_observations(as.integer(nsamples[s]), c(0), 
                                         censor=0.3, C_z=1)
            
            breslow <- fdp_breslow(list(obs_tree), c(epsilons[e]), c(0), p_hat, 
                                   C_z=1, delta=c(0.001), cutoff=1)
            logepscdp_NA[e, s, r] <-  sup_norm_diff(breslow$vals, breslow$times,
                                                    true_baseline=function(x){x}) 
        }
    }
}

# log plots for epsilon
line_plot(aperm(log(logepscdp_beta, 2), perm=c(2, 1, 3)),
          log(epsilons, 2), nsamples, bar=0.05, xlab=expression(log(epsilon)),
          colour='Samples', ylab='log(Squared error)')
ggsave('cdp_beta_logepsilon.png', width=4, height=3)

line_plot(aperm(log(logepscdp_breslow, 2), perm=c(2, 1, 3)),
          log(epsilons, 2), nsamples, bar=0.05, xlab=expression(log(epsilon)),
          colour='Samples', ylab='log(Sup error)')
ggsave('cdp_breslow_logepsilon.png', width=4, height=3)

line_plot(aperm(log(logepscdp_NA, 2), perm=c(2, 1, 3)),
          log(epsilons, 2), nsamples, bar=0.05, xlab=expression(log(epsilon)),
          colour='Samples', ylab='log(Sup error)')
ggsave('cdp_NA_logepsilon.png', width=4, height=3)

save(lognumcdp_beta, lognumcdp_breslow, lognumcdp_NA, logepscdp_beta, 
     logepscdp_breslow, logepscdp_NA, file='CDP_numsims.Rdata')
