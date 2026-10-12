library('FDPCox')
source('experiment_functions.R')

set.seed(123)

reps <- 100
epsilons <- c(0.2, 1, 5, 15)
nsamples <- 2^(12:18)
lognumcdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
lognumcdp_NA <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, sensitivity=NA)
        lognumcdp_beta[e, s, ] <- results$beta_error
        lognumcdp_NA[e, s, ] <- results$lambda_error
    }
}


epsilons <- 2^seq(-4, 10, 2)
nsamples <- c(10000, 20000, 50000)
logepscdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
logepscdp_NA <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, sensitivity=NA)
        logepscdp_beta[e, s, ] <- results$beta_error
        logepscdp_NA[e, s, ] <- results$lambda_error
    }
}



save(lognumcdp_beta, lognumcdp_NA, logepscdp_beta, logepscdp_NA, file='CDP_numsims.Rdata')

# log plots for samples
line_plot(log(lognumcdp_beta, 2), round(log(nsamples, 2), 2), epsilons, bar=0.05,
          ylab='log(Square error)', xlab='log(samples)')
ggsave('cdp_beta_logsamples.png', width=5, height=3)
line_plot(log(lognumcdp_NA, 2), round(log(nsamples, 2), 2), epsilons, bar=0.05,
          ylab='log(Sup error)', xlab='log(samples)')
ggsave('cdp_NA_logsamples.png', width=5, height=3)


# log plots for epsilon
line_plot(aperm(log(logepscdp_beta, 2), perm=c(2, 1, 3)),
          log(epsilons, 2), nsamples, bar=0.05, xlab=expression(log(epsilon)),
          colour='Samples', ylab='log(Squared error)')
ggsave('cdp_beta_logepsilon.png', width=5, height=3)

line_plot(aperm(log(logepscdp_NA, 2), perm=c(2, 1, 3)),
          log(epsilons, 2), nsamples, bar=0.05, xlab=expression(log(epsilon)),
          colour='Samples', ylab='log(Sup error)')
ggsave('cdp_NA_logepsilon.png', width=5, height=3)
