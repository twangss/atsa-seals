library(MARSS)
dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row

# n_states <- max(as.numeric(factor(custom_vec)))
# r_unequal <- seq(1, nrow(species_i_df))
# r_equal <- rep(1, nrow(species_i_df))
# uq_unequal <- seq(1, n_states)
# uq_equal <- rep(1, n_states)

f <- fit_stan(y = dat,
              mcmc_list = list(n_mcmc = 1000, n_burn = 500, n_chain = 1, n_thin = 1),
              marss = list(states = c(1,2,3), 
                           obsVariances = c(1,1,1), 
                           proVariances = c(1,2,3),
                           trends = c(1,2,3)))

#library(cmdstanr)
#file <- file.path("inst/stan/marss.stan")
#mod <- cmdstan_model(file)

