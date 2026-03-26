library('FDPCox')
source('experiment_functions.R')

set.seed(123)
reps <- 200
epsilons <- c(0.6, 0.8, 1, 1.5, 3, 6)
nservers <- 2*(1:6)

# experiments satisfying the FDP constraints 
fdps_beta <- array(0, dim=c(length(epsilons), length(nservers), reps))
fdps_hazard <- array(0, dim=c(length(epsilons), length(nservers), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nservers)){
        print(nservers[s])
        results <- fdp_experiment(S=nservers[s], nsamples=rep(20000, nservers[s]), 
                                  beta=c(0, 0.5, 0.8), 
                                  epsilon=rep(epsilons[e], nservers[s]), 
                                  censor=0.3, reps=reps, interactive=FALSE)
        fdps_beta[e, s, ] <- results$beta_error
        fdps_hazard[e, s, ] <- results$lambda_error
    }
}

# experiments where data is reused for fully-interactive mechanism. 
fdp_beta <- array(0, dim=c(length(epsilons), length(nservers), reps))
fdp_hazard <- array(0, dim=c(length(epsilons), length(nservers), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(nservers)){
        print(nservers[s])
        results <- fdp_experiment(nservers[s], rep(5000, nservers[s]), 
                                  c(0, 0.5, 0.8), rep(epsilons[e], nservers[s]), 
                                  censor=0.3, reps=reps, interactive=TRUE,
                                  trunc_mul=0.8)
        fdp_beta[e, s, ] <- results$beta_error
        fdp_hazard[e, s, ] <- results$lambda_error
    }
}

line_plot(fdps_beta, nservers, epsilons, bar=0.3, xlab='Servers')
ggsave('fdpalg_beta.png', width=5, height=3)
line_plot(fdps_hazard, nservers, epsilons, ylab='Sup error', xlab='Servers', bar=0.3)
ggsave('fdpalg_breslow.png', width=5, height=3)

line_plot(fdp_beta, nservers, epsilons, bar=0.3, xlab='Servers')
ggsave('fdpnt_beta.png', width=5, height=3)
line_plot(fdp_hazard, nservers, epsilons, bar=0.3, xlab='Servers', ylab='Sup error')
ggsave('fdpnt_breslow.png', width=5, height=3)


