#source('main_functions.R')
#source('helper_functions.R')
source('experiment_functions.R')

set.seed(123)
reps <- 200
epsilons <- c(0.6, 0.8, 1, 1.5, 3, 6)
censoring_rate <- seq(from=0.1, to=1.3, by=0.2)
censoring_beta <- array(0, dim=c(length(epsilons), length(censoring_rate), reps))
censoring_breslow <- array(0, dim=c(length(epsilons), length(censoring_rate), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(censoring_rate)){
        print(censoring_rate[s])
        results <- cdp_full_experiment(10000, c(0, 0.5, 0.8), epsilons[e], 
                                       reps=reps, censor=censoring_rate[s])
        censoring_beta[e, s, ] <- results$beta_error
        censoring_breslow[e, s, ] <- results$lambda_error
    }
}

save(censoring_beta, censoring_breslow, file='censoring.Rdata')

line_plot(censoring_beta, censoring_rate, epsilons, bar=0.05,
          xlab='Censoring distribution rate')
ggsave('censoring_beta.png', width=5, height=3)
line_plot(censoring_breslow, censoring_rate, epsilons, bar=0.05,
          xlab='Censoring distribution rate', ylab='Sup error')
ggsave('censoring_breslow.png', width=5, height=3)
