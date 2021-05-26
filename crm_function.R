

# Author: Dario Zocholl, Instiute of Biometry and Clinical Epidemiology, Charité - Universitätsmedizin Berlin

########################################
######## This script is required for the simulations for the manuscript "On the feasibility of pediatric dose-finding trials in small samples with information from a preceding trial in adults."
######## This script contains a function to simulate the conduct of phase I trials using the Continual Reassessment Method 
######## All files used for the manuscript are publicly available under https://github.com/dariozchl/Ped-PhaseI-Simulations
########################################


sim.phaseI <- function(algorithm, doses, target.tox, true.tox,
                       prior, borrow, prior.distr,
                       skeleton = NULL, dR = NULL, 
                       sample.size, n.sim, cohort.size, 
                       stopping.rule, escalation.rule,
                       iterations, stanmodel){ 
  
  start.time <- Sys.time() 
  
  library(rstan)
  rstan_options(auto_write = TRUE)
  
  # specify inv.logit-function()
  inv.logit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  
  # specify model and data input for the STAN model
  if (algorithm == "CRM" && borrow == "fixed" && prior.distr == "Normal"){
    
    mu <- prior[[1]]; sigma <- prior[[2]]
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$beta), ncol=1, byrow = F)
    }
    distr.tox.fun <- function(){skeleton[k]^exp(posterior)}
    mean.tox.fun <- function(){skeleton[k]^exp(mean(posterior))}
    
    
    } else if (algorithm == "CRM" && borrow == "partial" && prior.distr == "Normal"){
      
      mu <- prior[[1]]; sigma <- prior[[2]]
      data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton)}
      mcmc.sample <- function(){
        fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
        matrix(c(rstan::extract(fit)$beta), ncol=1, byrow = F)
      }
      distr.tox.fun <- function(){skeleton[k]^exp(posterior)}
      mean.tox.fun <- function(){skeleton[k]^exp(mean(posterior))}
      
    
  } else if (algorithm == "CRM" && borrow == "hierarchical" && prior.distr == "Normal"){
    
    mu <- prior[[1]]; sigma <- prior[[2]]
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$beta, rstan::extract(fit)$omega), ncol=2, byrow = F)
    }
    distr.tox.fun <- function(){skeleton[k]^exp(posterior[,1])}
    mean.tox.fun <- function(){skeleton[k]^exp(mean(posterior))}
    

    
  } else if (algorithm == "CRM" && borrow == "mixture" && prior.distr == "Normal"){

    mu1 <- prior[[1]][1]; sigma1 <- prior[[1]][2]; mu2 <- prior[[2]][1]; sigma2 <- prior[[2]][2]
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "mu1" = mu1, "mu2" = mu2, 
                                "sigma1" = sigma1, "sigma2" = sigma2, "skeleton" = skeleton)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$beta), ncol=1, byrow = F)
    }
    distr.tox.fun <- function(){skeleton[k]^exp(posterior)}
    mean.tox.fun <- function(){skeleton[k]^exp(mean(posterior))}
    
    
    
  } else if(algorithm == "PLM2" && borrow == "fixed") {
    
    mu <- c(prior[[1]],prior[[2]])
    Sigma <- prior[[3]]
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "dR"=dR, "mu"=mu, "Sigma"=Sigma)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$logalpha, rstan::extract(fit)$logbeta), ncol=2, byrow = F)
    }
    distr.tox.fun <- function(){inv.logit(posterior[,1] + exp(posterior[,2]) * log(doses[k]/dR))}
    mean.tox.fun <- function(){inv.logit(mean(posterior[,1]) + exp(mean(posterior[,2])) * log(doses[k]/dR))}
    
    
    
  } else if(algorithm == "PLM2" && borrow == "partial") {
    
    mu <- c(prior[[1]],prior[[2]])
    Sigma <- prior[[3]]
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "dR"=dR, "mu"=mu, "Sigma"=Sigma)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$logalpha, rstan::extract(fit)$logbeta), ncol=2, byrow = F)
    }
    distr.tox.fun <- function(){inv.logit(posterior[,1] + exp(posterior[,2]) * log(doses[k]/dR))}
    mean.tox.fun <- function(){inv.logit(mean(posterior[,1]) + exp(mean(posterior[,2])) * log(doses[k]/dR))}
    
    
    
  } else if(algorithm == "PLM2" && borrow == "hierarchical") {
    
    mu <- prior[[1]]; sigma <- prior[[2]]; rho <- prior[[3]]
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "dR"=dR, "mu"=mu, "sigma"=sigma, "rho" = rho)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$logalpha, rstan::extract(fit)$logbeta, rstan::extract(fit)$omega), ncol=3, byrow = F)
    }
    distr.tox.fun <- function(){inv.logit(posterior[,1] + exp(posterior[,2]) * log(doses[k]/dR))}
    mean.tox.fun <- function(){inv.logit(mean(posterior[,1]) + exp(mean(posterior[,2])) * log(doses[k]/dR))}

    
    
  } else if(algorithm == "PLM2" && borrow == "mixture") {

    mu_1 <- prior[[1]][[1]]; mu_2 <- prior[[2]][[1]]; Sigma_1 <- prior[[1]][[2]]; Sigma_2 <- prior[[2]][[2]]; 
    data.fun <- function(){list("K"=K, "tox"=tox, "N"=N, "doses"=doses, "dR"=dR, 
                                "mu_1"=mu_1, "mu_2"=mu_2, "Sigma_1"=Sigma_1, "Sigma_2"=Sigma_2)}
    mcmc.sample <- function(){
      fit <- sampling(stanmodel, data = data, warmup = 2000, iter = iterations, chains = 2, cores = 1, thin = 1, refresh = 0)
      matrix(c(rstan::extract(fit)$logalpha, rstan::extract(fit)$logbeta), ncol=2, byrow = F)
    }
    distr.tox.fun <- function(){inv.logit(posterior[,1] + exp(posterior[,2]) * log(doses[k]/dR))}
    mean.tox.fun <- function(){inv.logit(mean(posterior[,1]) + exp(mean(posterior[,2])) * log(doses[k]/dR))}
    
    
  } else {return("please specify the algorithm and type of borrowing correctly")}
  
  
  
  if (stopping.rule == "median") {
    stop.rule <- function(){median.tox[1] > 0.35}
  } else if(stopping.rule == "mean") {
    stop.rule <- function(){mean.tox[1] > (1/3)}
  } else {return("error")}
  
  if (escalation.rule == "median") {
    escalation.rule <- function(){ifelse(which.min(abs(median.tox - target.tox)) > k.star, k.star+1, which.min(abs(median.tox - target.tox)))} 
  } else if(escalation.rule == "mean") {
    escalation.rule <- function(){ifelse(which.min(abs(mean.tox - target.tox)) > k.star, k.star+1, which.min(abs(mean.tox - target.tox)))}
  } else {return("error")}
  
  
  
  
  
  results.list <- list()
  posterior.list <- list()
  
  
  for(nsim in 1:n.sim){ 
    
    K <- length(doses)
    tox <- c(rep(0,length(doses)))
    notox <- c(rep(0,length(doses)))
    tox.cohort <- c(rep(0, length(1:(sample.size/cohort.size))))
    mean.tox <- 0 # if mean.tox starts with only one entry, the first k.star will always be the first dose level 
    k.star <- 1
    distr.tox <- c()
    q95.tox <- c()
    q05.tox <- c()
    
    
    for(cohort in 1:(sample.size/cohort.size)){
      # treat cohort of patients and count toxicities
      tox.cohort[cohort] <- rbinom(1, cohort.size, true.tox[k.star])
      tox[k.star] <- tox[k.star] + tox.cohort[cohort]
      notox[k.star] <- notox[k.star] + cohort.size - tox.cohort[cohort]
      
      # count patients treated so far
      N <- tox + notox
      
      # sample model parameters from posterior distribution
      data <- data.fun()
      
      posterior <- mcmc.sample()
      
      
      # analysis of posterior toxicity probability for each dose
      for(k in 1:length(doses)){
        
        # construct posterior distribution of toxicity probability
        distr.tox[[k]] <- distr.tox.fun()
        
        # calculate mean toxicity
        mean.tox[k] <- mean.tox.fun()
      }
      
      # determine the dose for the following cohort
      k.star <- escalation.rule()
      
      # STOP criterion
      if (stop.rule()) k.star = 0
      if (stop.rule()) break
    }
    
    # define the MTD
    MTD <- k.star
    
    # # store weighting parameters (if weighting was applied)
    # omega <- ifelse(algorithm == "PLM2" && borrow == "hierarchical", mean(posterior[,3]), 
    #                 ifelse(borrow == "partial", omega, NA))
    # delta <- ifelse(algorithm == "CRM" && borrow == "hierarchical", mean(posterior[,2]), NA)
    
    
    
    # save the results
    if (nsim == 1) {
      results <- data.frame("dose.level" = 1:K, "doses" = doses, "mean.tox" = round(mean.tox,3), "true.tox" = true.tox, "toxicities.per.dose" = tox, "patients.per.dose" = N,
                            "MTD.dose.level" = k.star, "simulation.ID" = nsim, "sample.size" = sample.size)
    } else {
      results <- rbind(results, data.frame("dose.level" = 1:K, "doses" = doses, "mean.tox" = mean.tox, "true.tox" = true.tox, "toxicities.per.dose" = tox, "patients.per.dose" = N,
                                           "MTD.dose.level" = k.star, "simulation.ID" = nsim, "sample.size" = sample.size))             
    }
    
    if (nsim == 1) {
      posterior.df <- cbind(posterior, "simulation.ID"=1)
    } else {
      posterior.df <- rbind(posterior.df, cbind(posterior, "simulation.ID"=nsim))
    }
  }
  
  return(list(results, posterior.df))
}


