source('helper_functions.R')
source('main_functions.R')
source('experiment_functions.R')
library(survival)
set.seed(123)

rotterdam_covs <- rotterdam[c('hormon', 'chemo', 'meno')]
rotterdam_covs$grade <- as.integer(rotterdam$grade == 3)
rotterdam_covs$sizemed <- as.integer(rotterdam$size=='20-50')
rotterdam_covs$sizelarge <- as.integer(rotterdam$size=='>50')

rotterdam_covs <- rotterdam_covs/sqrt(5)
rotterdam_covs <- as.matrix(rotterdam_covs)
rfs  <- pmax(rotterdam$recur, rotterdam$death)
rfstime <- with(rotterdam, ifelse(recur==1, rtime, dtime))
non_priv_fit <- coxph(Surv(rfstime, rfs) ~.,data=data.frame(rotterdam_covs))
non_priv_coefs <- non_priv_fit$coefficients
#Breslow at covariates=0
non_priv_breslow <- basehaz(non_priv_fit, centered=FALSE)

print(non_priv_coefs)
print(sqrt(sum(non_priv_coefs^2)))

epsilons <- c(3, 6, 9, 12, 15)
reps <- 10
total_samples <- nrow(rotterdam_covs)
servers <- 3
n <- total_samples%/%servers
real_coefs <- array(0, dim=c(servers, length(epsilons), reps))
real_breslow <- array(0, dim=c(servers, length(epsilons), reps))
real_semiprivbreslow <- array(0, dim=c(servers, length(epsilons), reps))


for (rep in 1:reps){
    if (rep%%40==0){print(rep)}
    # construct server datasets
    mixed_indices <- sample(nrow(rotterdam_covs))
    for (numservers in 1:servers){
        obs_list <- list()
        for (s in 1:numservers){
        server_indices <- mixed_indices[((s-1)*n+1):(s*n)]
        obs_list[[s]] <- cbind(rfstime[server_indices], rfs[server_indices], 
                               rotterdam_covs[server_indices, ])
        rownames(obs_list[[s]]) <- NULL
        }
        
        for (e in 1:length(epsilons)){
        epsilon <- rep(epsilons[e], numservers)
        # Cox regression
        priv_coefs <- fdp_cox_interactive(obs_list, epsilon=epsilon,
                                          C_beta=3, cutoff=7000, stepsize=0.5,
                                          niters=300)
        real_coefs[numservers, e, rep] <- sum((priv_coefs - non_priv_coefs)^2)/
                                                sum(non_priv_coefs^2)
        # Breslow
        times_list <- list()
        for(s in 1:numservers){
            times_list[[s]] <- obs_list[[s]][, 1]
        }
        p_hat <- fdp_probabilities(times_list, epsilon, cutoff=3500)
        priv_breslow <- fdp_breslow(obs_list, epsilon, priv_coefs, p_hat, cutoff=3500,
                                    truncation=p_hat*0.9)
        semipriv_breslow <- fdp_breslow(obs_list, epsilon, non_priv_coefs, 
                                           p_hat, cutoff=3500, truncation=p_hat*0.9)
        real_breslow[numservers, e, rep] <- step_functions_sup_diff(priv_breslow$times, 
                                                              priv_breslow$vals,
                                                              non_priv_breslow$time,
                                                              non_priv_breslow$hazard)
        real_semiprivbreslow[numservers, e, rep] <- step_functions_sup_diff(semipriv_breslow$times, 
                                                              semipriv_breslow$vals,
                                                              non_priv_breslow$time,
                                                              non_priv_breslow$hazard)
        }
    }
}

save.image(file='real_experiments.RData')

real_CDPcoefs <- matrix(0, length(epsilons), reps)
real_CDPbreslow <- matrix(0, length(epsilons), reps)
real_CDPsemipriv <- matrix(0, length(epsilons), reps)
obs_list <- cbind(rfstime, rfs, rotterdam_covs)
for (rep in 1:reps){
    if (rep%%25==0){print(rep)}
    for (e in 1:length(epsilons)){
        priv_coefs <- cdp_cox(obs_list, epsilon=epsilons[e], C_beta=3, cutoff=7000, 
                              stepsize=0.5, niters=300)
        real_CDPcoefs[e, rep] <- sum((priv_coefs - non_priv_coefs)^2)/
            sum(non_priv_coefs^2)
        times_list <- list(obs_list[, 1])
        p_hat <- fdp_probabilities(times_list, epsilon=epsilons[e], cutoff=3500)
        priv_breslow <- fdp_breslow(list(obs_list), epsilon=epsilons[e], priv_coefs, 
                                    p_hat, cutoff=3500, truncation=p_hat*0.9)
        semipriv_breslow <- fdp_breslow(list(obs_list), epsilons[e], non_priv_coefs, 
                                        p_hat, cutoff=3500, truncation=p_hat*0.9)
        real_CDPbreslow[e, rep] <- step_functions_sup_diff(priv_breslow$times, 
                                                        priv_breslow$vals,
                                                        non_priv_breslow$time,
                                                        non_priv_breslow$hazard)
        real_CDPsemipriv[e, rep] <- step_functions_sup_diff(semipriv_breslow$times, 
                                                            semipriv_breslow$vals,
                                                            non_priv_breslow$time,
                                                            non_priv_breslow$hazard)
    }
}

real_CDPcoefs <- array(real_CDPcoefs, dim=c(1, dim(real_CDPcoefs)))
real_CDPbreslow <- array(real_CDPbreslow, dim=c(1, dim(real_CDPbreslow)))
real_CDPsemipriv <- array(real_CDPsemipriv, dim=c(1, dim(real_CDPsemipriv)))
real_coefs <- abind::abind(real_coefs, real_CDPcoefs, along=1)
real_breslow <- abind::abind(real_breslow, real_CDPbreslow, along=1)
real_semiprivbreslow <- abind::abind(real_semiprivbreslow, real_CDPsemipriv, along=1)

line_plot(real_coefs, epsilons, c(1, 2, 3, 'CDP'), bar=0.5, 
          xlab=expression(epsilon), colour='Servers',
          ylab='Standardised squared error')
ggsave('real_coefs.png', width=5, height=4)
line_plot(real_breslow, epsilons, c(1, 2, 3, 'CDP'), bar=0.5, 
          xlab=expression(epsilon), colour='Servers', ylab='Sup error')
ggsave('real_breslow.png', width=5, height=4)
line_plot(real_semiprivbreslow, epsilons, c(1, 2, 3, 'CDP'), bar=0.5, 
          xlab=expression(epsilon), colour='Servers', ylab='Sup error')
ggsave('real_semiprivbreslow.png', width=5, height=4)



save(real_CDPcoefs, real_CDPbreslow, real_CDPsemipriv, file='realCDP.RData')
save(real_coefs, real_breslow, real_semiprivbreslow, file='real.RData')

