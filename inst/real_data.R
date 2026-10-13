library('FDPCox')
source('experiment_functions.R')
library(survival)
library(abind)
set.seed(123)

rotterdam_covs <- rotterdam[c('hormon', 'chemo', 'meno')]
rotterdam_covs$grade <- as.integer(rotterdam$grade == 3)
rotterdam_covs$sizemed <- as.integer(rotterdam$size=='20-50')
rotterdam_covs$sizelarge <- as.integer(rotterdam$size=='>50')

rotterdam_covs <- rotterdam_covs/sqrt(5)
rotterdam_covs <- as.matrix(rotterdam_covs)
rfs  <- pmax(rotterdam$recur, rotterdam$death)
rfstime <- with(rotterdam, ifelse(recur==1, rtime, dtime))
non_priv_fit <- coxph(
    Surv(pmin(rfstime, 7000),
         as.integer(rfs == 1 & rfstime <= 7000)) ~ .,
    data = data.frame(rotterdam_covs),
    ties = "breslow"
)
non_priv_coefs <- non_priv_fit$coefficients
# Breslow estimator at covariates=0
non_priv_breslow <- basehaz(non_priv_fit, centered=FALSE)
print(non_priv_coefs)
print(sqrt(sum(non_priv_coefs^2)))

na_fit <- survfit(Surv(rfstime, rfs) ~ 1, ctype = 1)
na_est <- data.frame(time = na_fit$time, hazard = na_fit$cumhaz)
kept_times <- na_est$time>0
C_lambda_hat <- 3*max(na_est$hazard[kept_times]/na_est$time[kept_times])

epsilons <- c(3, 6, 9, 12, 15)
reps <- 200
total_samples <- nrow(rotterdam_covs)
servers <- 3
n <- total_samples%/%servers
real_coefs <- array(0, dim=c(servers, length(epsilons), reps))
real_NA <- array(0, dim=c(servers, length(epsilons), reps))


for (rep in 1:reps){
    # construct server datasets
    mixed_indices <- sample(nrow(rotterdam_covs))
    for (numservers in 1:servers){
        obs_list_beta <- list()
        obs_list_na <- list()
        for (s in 1:numservers){
            server_indices <- mixed_indices[((s-1)*n+1):(s*n)]
            obs_list_beta[[s]] <- cbind(rfstime[server_indices], rfs[server_indices], 
                                        rotterdam_covs[server_indices, ])
            obs_list_na[[s]] <- cbind(rfstime[server_indices], rfs[server_indices])
        }
        
        for (e in 1:length(epsilons)){
            epsilon <- rep(epsilons[e], numservers)
            # Cox regression
            priv_coefs <- fdp_cox_interactive(obs_list_beta, epsilon=epsilon, C_beta=3, 
                                              cutoff=7000, stepsize=0.5, niters=200)
            real_coefs[numservers, e, rep] <- sum((priv_coefs - non_priv_coefs)^2)/
                sum(non_priv_coefs^2)
            # NA
            priv_na <- fdp_breslow(obs_list_na, epsilon, cutoff=3500, C_lambda=C_lambda_hat)
            real_NA[numservers, e, rep] <- step_functions_sup_diff(priv_na$times,priv_na$vals,
                                                                   na_est$time, na_est$hazard,
                                                                   cutoff=3500)
        }
    }
}

# run the CDP versions of the experiment, using full data
real_CDPcoefs <- matrix(0, length(epsilons), reps)
real_CDP_na <- matrix(0, length(epsilons), reps)
obs_list_beta <- cbind(rfstime, rfs, rotterdam_covs)
obs_list_na <- cbind(rfstime, rfs)

for (rep in 1:reps){
    obs_na_rep <- obs_list_na[sample.int(nrow(obs_list_na)), , drop=FALSE]
    for (e in 1:length(epsilons)){
        priv_coefs <- cdp_cox(obs_list_beta, epsilon=epsilons[e], C_beta=3, cutoff=7000, 
                              stepsize=0.5, niters=200)
        real_CDPcoefs[e, rep] <- sum((priv_coefs - non_priv_coefs)^2)/
            sum(non_priv_coefs^2)
        priv_na <- fdp_breslow(list(obs_na_rep), epsilon=c(epsilons[e]), cutoff=3500, 
                               C_lambda=C_lambda_hat)
        real_CDP_na[e, rep] <- step_functions_sup_diff(priv_na$times, priv_na$vals,
                                                       na_est$time, na_est$hazard,
                                                       cutoff=3500)
    }
}

save(real_coefs, real_NA, real_CDPcoefs, real_CDP_na, file='real_experiments.RData')

real_CDPcoefs <- array(real_CDPcoefs, dim=c(1, dim(real_CDPcoefs)))
real_CDP_na <- array(real_CDP_na, dim=c(1, dim(real_CDP_na)))
real_coefs <- abind::abind(real_coefs, real_CDPcoefs, along=1)
real_NA <- abind::abind(real_NA, real_CDP_na, along=1)

line_plot(real_coefs, epsilons, c(1, 2, 3, 'CDP'), bar=0.5, 
          xlab=expression(epsilon), colour='Servers',
          ylab='Standardised squared error')
ggsave('real_coefs.png', width=5, height=3)
line_plot(real_NA, epsilons, c(1, 2, 3, 'CDP'), bar=0.5, 
          xlab=expression(epsilon), colour='Servers', ylab='Sup error')
ggsave('real_NA.png', width=5, height=3)
