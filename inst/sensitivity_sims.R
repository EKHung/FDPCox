#source('main_functions.R')
#source('helper_functions.R')
source('experiment_functions.R')

set.seed(123)
reps <- 200
epsilons <- c(0.6, 0.8, 1, 1.5, 3, 6)

Cz_values <- seq(from=0.8, to=2, by=0.2)
sens_sensitivity <- array(0, dim=c(length(epsilons), length(Cz_values), reps))
for (e in 1:length(epsilons)){
    for (s in 1:length(Cz_values)){
        sens_sensitivity[e, s, ] <- cdp_beta_experiment(10000, c(0, 0.5, 0.8),
                                                        epsilon=epsilons[e], 
                                                        reps=reps, 
                                                        sens_Cz=Cz_values[s])
    }
}

step <- seq(from=0.2, to=0.8, by=0.1)
step_sensitivity <- array(0, dim=c(length(epsilons), length(step), reps))
for (e in 1:length(epsilons)){
    print(epsilons[e])
    for (s in 1:length(step)){
        step_sensitivity[e, s, ] <- cdp_beta_experiment(10000, c(0, 0.5, 0.8),
                                                        epsilon=epsilons[e], 
                                                        reps=reps,
                                                        stepsize=step[s])
    }
}

save(sens_sensitivity, step_sensitivity, file='sensitivity.Rdata')

line_plot(sens_sensitivity, Cz_values, epsilons, bar=0.05, 
          xlab='Covariate bound')
ggsave('noise_sensitivity.png', width=5, height=3)

line_plot(step_sensitivity, step, epsilons, bar=0.03, xlab='Step size')
ggsave('stepsize_sensitivity.png', width=5, height=3)


