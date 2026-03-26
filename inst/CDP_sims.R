library('FDPCox')
source('experiment_functions.R')

set.seed(123)
epsilons <- c(0.6, 0.8, 1, 1.5, 3, 6)
nsamples <- seq(from=2000, to=14000, by=2000)
reps <- 200
labelcdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
labelcdp_breslow <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, sensitivity='label')
        labelcdp_beta[e, s, ] <- results$beta_error
        labelcdp_breslow[e, s, ] <- results$lambda_error
    }
}

cdp_beta <- array(0, dim=c(length(epsilons), length(nsamples), reps))
cdp_breslow <- array(0, dim=c(length(epsilons), length(nsamples), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nsamples)){
        results <- cdp_full_experiment(nsamples[s], c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, sensitivity=NA)
        cdp_beta[e, s, ] <- results$beta_error
        cdp_breslow[e, s, ] <- results$lambda_error
    }
}

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
            obs <- sim_observations(10000, beta, censor=0.3)
            out <- cdp_cox(obs, epsilons[s])
            cdp_dimension[e, s, r] <- sum((beta - out)^2)
        }
    }
}
cdp_dimension <- aperm(cdp_dimension, perm = c(2, 1, 3))

save(labelcdp_beta, labelcdp_breslow, cdp_beta, cdp_breslow, cdp_dimensions, 
     file='CDP_sims.Rdata')


line_plot(cdp_beta, nsamples, epsilons, bar=400)
ggsave('cdp_beta.png', width=4.5, height=3)
line_plot(cdp_breslow, nsamples, epsilons, ylab='Sup error', bar=400)
ggsave('cdp_breslow.png', width=4.5, height=3)
line_plot(cdp_dimension, param_dim, epsilons, xlab='Dimension', bar=0.2)
ggsave('cdp_dimension.png', width=4.5, height=3)



labelcdp_beta_df <- line_plot(labelcdp_beta, nsamples, epsilons, df_func=T)
cdp_beta_df <- line_plot(cdp_beta, nsamples, epsilons, df_func=T)

ggplot() +
    geom_line(data=labelcdp_beta_df,
              mapping=aes(x=x, y=mean, colour=s, linetype = "Label")) +
    geom_errorbar(data=labelcdp_beta_df,
                  mapping=aes(x=x, ymin=lower, ymax=upper, colour=s,
                              linetype = "Label", width=400)) +
    geom_line(data=cdp_beta_df,
              mapping=aes(x=x, y=mean, colour=s, linetype = "Full")) +
    geom_errorbar(data=cdp_beta_df,
                  mapping=aes(x=x, ymin=lower, ymax=upper, colour=s,
                              linetype = "Full", width=550)) +
    scale_linetype_manual(values=c(Full="dashed", Label="solid"), name="") +
    scale_x_continuous(breaks=nsamples) +
    labs(x='Samples', y='Squared error', colour=expression(epsilon)) +
    theme_minimal()
ggsave('labelcomp_beta.png', width=5, height=4)

