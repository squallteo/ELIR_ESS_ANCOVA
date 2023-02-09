rm(list=ls())
library(RBesT)
library(rjags)
library(doParallel)
library(mvtnorm)

#design parameters
n_total <- 120
RRatio <- 2 #randomization ratio: treatment:control
ped_means_pbo <- c(8.0, 8.0)
ped_means_active <- c(8.0, 8.0)
ped_vcov <- matrix(c(1.0^2, 0.4, 0.4, 1.43^2), nrow=2)

prior_m <- -0.6
prior_s <- 0.6
prior_mix <- mixnorm(c(1, prior_m, prior_s), param="ms")

jags.script <- "
  model{
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1/ped_var)
      mu[i] <- beta0 + beta1*x[i,1] + delta_p*x[i,3]
    }
    
    beta0 ~ dnorm(0, 1/10000)
    beta1 ~ dnorm(0, 1/10000)
    
    delta_p ~ dnorm(adult_mean, 1/adult_var)
  }
"

ncores <- min(parallel::detectCores(), 40)
cl <- makeCluster(ncores-1)
registerDoParallel(cl)
#############################################
#############################################
#############################################
#parallel computing
results <- 
  
  foreach(s = 1:500, .combine = rbind, .packages = c("RBesT", "mvtnorm", "rjags"), .errorhandling = "remove") %dopar% {
    
    set.seed(s+712)
    
    pbo_dat <- rmvnorm((1/(1+RRatio))*n_total, ped_means_pbo, ped_vcov)
    active_dat <- rmvnorm((RRatio/(1+RRatio))*n_total, ped_means_active, ped_vcov)
    total_dat <- rbind(cbind(pbo_dat, rep(0, (1/(1+RRatio))*n_total)), cbind(active_dat, rep(1, (RRatio/(1+RRatio))*n_total)))
    alldt <- as.data.frame(cbind(total_dat, total_dat[,2]-total_dat[,1]))
    colnames(alldt) <- c('Base', "Aval", 'TRTP', 'Chg')
    
    #prior ESS of delta_p
    sigma(prior_mix) <- sqrt(2)*sd(alldt$Chg)
    prior_ess <- c(ess(prior_mix, "elir"),
                   ess(prior_mix, "morita"),
                   ess(prior_mix, "moment")
    )
    
    #run the Bayesian ANCOVA
    jags_data <- list(x = alldt[,1:3], y = alldt[,4], 
                      n = nrow(alldt),
                      adult_mean = prior_m,
                      adult_var = prior_s^2,
                      ped_var = var(alldt$Chg)
    )
    
    jag_obj <- jags.model(textConnection(jags.script),
                          data = jags_data,
                          n.chains = 3, 
                          n.adapt = 2000
    )
    
    jags_sample <- coda.samples(model=jag_obj,
                                variable.names = c('delta_p'),
                                n.iter=12000, 
                                thin = 2
    )
    
    tt<- unlist(jags_sample)
    
    #posterior ESS of delta_p
    post_hat <- automixfit(tt, type="norm")
    sigma(post_hat) <- sqrt(2)*sd(alldt$Chg)
    post_ess <- c(ess(post_hat, "elir"),
                  ess(post_hat, "morita"),
                  ess(post_hat, "moment")
    )
    
    c(prior_ess, post_ess)
    
  }

results <- as_tibble(results)
names(results) <- c("ELIR_prior", "Morita_prior", "Moment_prior",
                    "ELIR_post", "Morita_post", "Moment_post")
tt <- 
  results %>% mutate(ELIR_diff = ELIR_post - ELIR_prior,
                     Morita_diff = Morita_post - Morita_prior,
                     Moment_diff = Moment_post - Moment_prior)

#check the average differences between posterior and prior ESS
summary(tt$ELIR_diff)
summary(tt$Morita_diff)
summary(tt$Moment_diff)

# lm_obj <- summary(lm(Chg ~ Base + TRTP, lm_dat))
# rse <- lm_obj$sigma
