
# Author: Dario Zocholl, Instiute of Biometry and Clinical Epidemiology, Charité - Universitätsmedizin Berlin

########################################
######## This script performs the simulations for the manuscript "On the feasibility of pediatric dose-finding trials in small samples with information from a preceding trial in adults."
######## This script uses a function defined in the file "crm_function.R"
######## This script only simulates data - the analysis of the data is performed within the file "plots.R"
######## All files are publicly available under https://github.com/dariozchl/Ped-PhaseI-Simulations
########################################


library(rstan)
library(doParallel)
library(dplyr)
library(tidyr)
library(tibble)

set.seed(16)

setwd("S:/C01/iBikE/Studien/Phase1Studien/3_Programme/Data")
source("S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/crm_function.R")


########################################
######## Define the STAN-Models
########################################

generate_stan_model <- function(algorithm, borrow){
  # specify model and data input for the STAN model
  if (algorithm == "CRM" && borrow == "fixed"){
    smodel = "data {
      int<lower=1> K; // total number of dose levels 
      int<lower=0> N[K]; // total number of observations (so far) per dose level
      int<lower=0> tox[K]; // number of tox per dose
      vector[K] skeleton;
      real mu;
      real<lower=0> sigma;
    }
      
      
      parameters {
      real beta;
      }
      
      transformed parameters {
      real<lower=0,upper=1> p[K];
      
      for(k in 1:K){
      p[k] = pow(skeleton[k],exp(beta));
      }
      
      }
      
      model {
      tox ~ binomial(N,p);
      beta ~ normal(mu, sigma);
      }
      "
    return(stan_model(model_code = smodel))
    
  } else if (algorithm == "CRM" && borrow == "partial"){
    smodel = "data {
      int<lower=1> K; // total number of dose levels 
      int<lower=0> N[K]; // total number of observations (so far) per dose level
      int<lower=0> tox[K]; // number of tox per dose
      vector[K] skeleton;
      real mu;
      real<lower=0> sigma;
    }
      
      
      parameters {
      real beta;
      }
      
      transformed parameters {
      real<lower=0,upper=1> p[K];
      
      for(k in 1:K){
      p[k] = pow(skeleton[k],exp(beta));
      }
      
      }
      
      model {
      tox ~ binomial(N,p);
      beta ~ normal(mu, sigma);
      }
      "
    return(stan_model(model_code = smodel))
    
    
  } else if (algorithm == "CRM" && borrow == "hierarchical"){
    smodel = "data {
        int<lower=1> K; // total number of dose levels 
        int<lower=0> N[K]; // total number of observations (so far) per dose level
        int<lower=0> tox[K]; // number of tox per dose
        vector[K] skeleton;
        real mu;
        real<lower=0> sigma;
        }
        
        
        parameters {
        real beta;
        real<lower=0,upper=1> omega;
        }
        
        transformed parameters {
        real<lower=0,upper=1> p[K];
        real<lower=0> weighted_sigma; 

        for(k in 1:K){
        p[k] = pow(skeleton[k],exp(beta));
        }

        weighted_sigma = sigma/omega;
        
        }
        
        model {
        tox ~ binomial(N,p);
        beta ~ normal(mu, weighted_sigma);
        omega ~ beta(0.5,0.5);
        }
        "
    return(stan_model(model_code = smodel))
    
  } else if (algorithm == "CRM" && borrow == "mixture"){
    smodel = "data {
      int<lower=1> K; // total number of dose levels 
      int<lower=0> N[K]; // total number of observations (so far) per dose level
      int<lower=0> tox[K]; // number of tox per dose
      vector[K] skeleton;
      real mu1;
      real mu2;
      real<lower=0> sigma1;
      real<lower=0> sigma2;
      }
      
      
      parameters {
      real beta;
      }
      
      transformed parameters {    
        real<lower=0,upper=1> p[K];
        for (k in 1:K) {
          p[k] = pow(skeleton[k],exp(beta));
        }
      }
      
      model {
      target += log_sum_exp(normal_lpdf(beta | mu1, sigma1),
                            normal_lpdf(beta | mu2, sigma2)); 
      for (k in 1:K) {
      target += binomial_lpmf(tox[k] | N[k], pow(skeleton[k],exp(beta)));      
      }
      }
      "
    return(stan_model(model_code = smodel))
    
    
  } else if(algorithm == "PLM2" && borrow == "fixed") {
    smodel = "data {
        int<lower=1> K; // total number of dose levels 
        vector<lower=0>[K] doses;
        int<lower=0> tox[K]; // number of tox per dose
        int<lower=0> N[K]; // total number of patients per dose level
        int<lower=0> dR; // reference dose
        vector[2] mu;
        cov_matrix[2] Sigma;
        }
        
        parameters {          
          vector[2] logalphabeta;
        }
        
        transformed parameters {
    	    real<lower=0,upper=1> p[K];
          real logalpha;
          real logbeta;

          logalpha = logalphabeta[1];
          logbeta = logalphabeta[2];
          
          for(k in 1:K){
            p[k] = inv_logit(logalpha + exp(logbeta) * log(doses[k]/dR));
          }
        }
        
          
        model {
          logalphabeta ~ multi_normal(mu, Sigma);
          tox ~ binomial_logit(N, logalpha + exp(logbeta)*log(doses/dR));
        }
        "
    return(stan_model(model_code = smodel))
    
  } else if(algorithm == "PLM2" && borrow == "partial") {
    smodel = "data {
        int<lower=1> K; // total number of dose levels 
        vector<lower=0>[K] doses;
        int<lower=0> tox[K]; // number of tox per dose
        int<lower=0> N[K]; // total number of patients per dose level
        int<lower=0> dR; // reference dose
        vector[2] mu;
        cov_matrix[2] Sigma;
        }
        
        parameters {          
          vector[2] logalphabeta;
        }
        
        transformed parameters {
          real logalpha;
          real logbeta;
          
          logalpha = logalphabeta[1];
          logbeta = logalphabeta[2];
        }
          
        model {
          logalphabeta ~ multi_normal(mu, Sigma);
          tox ~ binomial_logit(N, logalpha + exp(logbeta)*log(doses/dR));
        }
        "
    return(stan_model(model_code = smodel))
    
    
  } else if(algorithm == "PLM2" && borrow == "hierarchical") {
    smodel = "data {
          int<lower=1> K; // total number of dose levels 
          vector<lower=0>[K] doses; // predictor variable 
          int<lower=0> tox[K]; // number of tox per dose
          int<lower=0> N[K]; // total number of patients per dose level
          int<lower=0> dR; // reference dose
          vector[2] mu;
          vector<lower=0>[2] sigma; 
          real rho;
          }
          
          parameters {
          vector[2] logalphabeta;
          real<lower=0,upper=1> omega;
          }
          
          transformed parameters {
          real logalpha;
          real logbeta;
          cov_matrix[2] Sigma;
          
          Sigma[1,1] = square(sigma[1])/omega;
          Sigma[2,2] = square(sigma[2])/omega;
          Sigma[1,2] = sigma[1]*sigma[2]*rho;
          Sigma[2,1] = sigma[1]*sigma[2]*rho;
          
          logalpha = logalphabeta[1];
          logbeta = logalphabeta[2];
          }
          
          model {
          logalphabeta ~ multi_normal(mu, Sigma);
          tox ~ binomial_logit(N, logalpha + exp(logbeta)*log(doses/dR));
          omega ~ beta(0.5,0.5);
          }
          "
    return(stan_model(model_code = smodel))
    
    
    
  } else if(algorithm == "PLM2" && borrow == "mixture") {
    smodel = "data {
          int<lower=1> K; // total number of dose levels 
          vector<lower=0>[K] doses; // predictor variable 
          int<lower=0> tox[K]; // number of tox per dose
          int<lower=0> N[K]; // total number of patients per dose level
          int<lower=0> dR; // reference dose
          vector[2] mu_1;
          vector[2] mu_2;
          cov_matrix[2] Sigma_1;
          cov_matrix[2] Sigma_2;
          }
          
          parameters {
          vector[2] logalphabeta;
          }

    	    transformed parameters {
    	    real<lower=0,upper=1> p[K];
          real logalpha;
          real logbeta;

          logalpha = logalphabeta[1];
          logbeta = logalphabeta[2];
          
          for(k in 1:K){
            p[k] = inv_logit(logalpha + exp(logbeta) * log(doses[k]/dR));
            }
          }

          model {
          target += log_sum_exp(multi_normal_lpdf(logalphabeta | mu_1, Sigma_1),
                                multi_normal_lpdf(logalphabeta | mu_2, Sigma_2)); 
          for (k in 1:K) {
            target += binomial_logit_lpmf(tox[k] | N[k], logalpha + exp(logbeta)*log(doses[k]/dR));
          }
          }
          "
    return(stan_model(model_code = smodel))
    
  } else {return("please specify the algorithm and type of borrowing correctly")}
  
}

crm.fixed <- generate_stan_model(algorithm = "CRM", borrow = "fixed")
crm.mixture <- generate_stan_model(algorithm = "CRM", borrow = "mixture")
#crm.hierarchical <- generate_stan_model(algorithm = "CRM", borrow = "hierarchical") # not used in the simulations
PLM2.fixed <- generate_stan_model(algorithm = "PLM2", borrow = "fixed")
PLM2.mixture <- generate_stan_model(algorithm = "PLM2", borrow = "mixture")
#PLM2.partial <- generate_stan_model(algorithm = "PLM2", borrow = "partial") # not used in the simulations



########################################
######## Define settings for adult trials
########################################



# define inv.logit-function()
inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

# define logit-function()
logit <- function(p){ 
  log(p / (1-p)) 
}

# define sample size for adults
ss.adults <- 30

# define doses
doses.adults <- c(5,10,15,20,27,36,47)


################################
# toxicity scenarios

LL4 <- function(x,b,c,d,e){c + (d-c)/(1+exp(b*(log(x)-log(e))))}
x=seq(1,100,by=0.1)

dose.tox.weak <- function(doses){1-LL4(x=doses,b=1.5,c=0,d=1,e=70)}
dose.tox.moderate <- function(doses){1-LL4(x=doses,b=1.3,c=0,d=1,e=50)}
dose.tox.strong <- function(doses){1-LL4(x=doses,b=1.8,c=0,d=1,e=27)}

probs.plateau <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)
spline.plateau <- splinefun(c(3.5, doses.adults,188,40), probs.plateau, "monoH.FC")

probs.waves <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)
spline.waves <- splinefun(c(doses.adults,188), probs.waves, "monoH.FC")

tox.adults.weak <- dose.tox.weak(doses.adults)
tox.adults.moderate <- dose.tox.moderate(doses.adults)
tox.adults.strong <- dose.tox.strong(doses.adults)
tox.adults.spline.plateau <- spline.plateau(doses.adults)
tox.adults.spline.waves <- spline.waves(doses.adults)

tox.adults <- list(tox.adults.weak, tox.adults.moderate, tox.adults.strong, tox.adults.spline.plateau, tox.adults.spline.waves)

# show all adult scenarios
round.2 <- function(x){round(x, 2)}
lapply(tox.adults, round.2)


# pediatric toxicity scenarios (t is an index for the toxicity scenario)
fun.true.tox.ped.weaker <- function(doses,t){ if(t==1){dose.tox.weak(doses)^1.3} else if(t==2){dose.tox.moderate(doses)^1.3} else if(t==3){dose.tox.strong(doses)^1.3} else if(t==4){spline.plateau.weaker(doses)} else if(t==5){spline.waves.weaker(doses)} }
fun.true.tox.ped.same <- function(doses,t){ if(t==1){dose.tox.weak(doses)} else if(t==2){dose.tox.moderate(doses)} else if(t==3){dose.tox.strong(doses)} else if(t==4){spline.plateau(doses)} else if(t==5){spline.waves(doses)} }
fun.true.tox.ped.stronger <- function(doses,t){ if(t==1){dose.tox.weak(doses)^0.7} else if(t==2){dose.tox.moderate(doses)^0.7} else if(t==3){dose.tox.strong(doses)^0.4} else if(t==4){spline.plateau.stronger(doses)} else if(t==5){spline.waves.stronger(doses)} }
spline.plateau.stronger <- function(doses){spline.plateau(doses)^0.7}
spline.plateau.weaker <- function(doses){spline.plateau(doses)^1.3}
spline.waves.stronger <- function(doses){spline.waves(doses)^0.7}
spline.waves.weaker <- function(doses){spline.waves(doses)^1.3}


################### Adult Prior ###################

# reference dose
dR.adults <- 20

# specify weakly informative prior for 2P-CRM
prior.PLM2 <- list(mu=c(logit(0.25), 0.4606), Sigma=matrix(c(1.939^2, 0, 0, 1.5^2), nrow = 2, byrow = T))

# show the priors for 1P- and 2P-CRM
# Pr(pi <= 0.8) = 0.9 for reference dose in 2P-CRM
inv.logit(qnorm(p=0.9,mean=-1.098612,sd=1.939))
# Pr(pi <= 0.8) = 0.9 for reference dose in 1P-CRM
0.25^exp(qnorm(p=0.1,mean=0,sd=1.425))

# expected odds ratio when dose is halved
exp(logit(inv.logit(logit(0.25) + exp(qnorm(p=0.5, mean=0.4605607, sd=1.5)) * log(10/20))) - logit(0.25))
# 10%-quantile of odds ratio when dose is halved
exp(logit(inv.logit(logit(0.25) + exp(qnorm(p=0.9, mean=0.4605607, sd=1.5)) * log(10/20))) - logit(0.25))
# 90%-quantile of odds ratio when dose is halved
exp(logit(inv.logit(logit(0.25) + exp(qnorm(p=0.1, mean=0.4605607, sd=1.5)) * log(10/20))) - logit(0.25))



########################################
######## Simulation adult trials
########################################

n.sim <- 1000
MTDs.adults <- list()
data.adults <- list()
start.time <- Sys.time()
posterior <- list()
output <- list(list(), list())


for(a in 1:length(tox.adults)){    
  start.time.local <- Sys.time()
  index <- a
  
  for(simulation.ID in 1:n.sim){
    output.temporary <- NULL
    while(!any(output.temporary[[1]]$mean.tox > (1/6) & output.temporary[[1]]$mean.tox < (1/3))){
      output.temporary <- sim.phaseI(algorithm="PLM2", doses = doses.adults, target.tox = 0.25, 
                                     true.tox = tox.adults[[a]],
                                     prior = list(unlist(prior.PLM2$mu)[1], unlist(prior.PLM2$mu)[2], prior.PLM2$Sigma), 
                                     borrow = "fixed", prior.distr = "Normal",
                                     skeleton = NULL, dR = dR.adults, sample.size = ss.adults, cohort.size = 3,
                                     n.sim = 1, iterations = 10000, stanmodel = PLM2.fixed,
                                     stopping.rule = "mean", escalation.rule = "mean")
    }  
    output.temporary[[2]] <- as_tibble(output.temporary[[2]])
    output.temporary[[1]]$simulation.ID <- output.temporary[[2]]$simulation.ID <- simulation.ID
    output[[1]] <- rbind(output[[1]], output.temporary[[1]] %>% add_column("tox.scenario"=a))
    output[[2]] <- rbind(output[[2]], output.temporary[[2]] %>% add_column("tox.scenario"=a))
  }
  
  cat("\nThis was adult trial no.", index, 
      ". Run time for this algorithm was ", round(difftime(Sys.time(), start.time.local, units = "mins"), 2), " minutes.",
      "\nTotal run time so far is ", round(difftime(Sys.time(), start.time, units = "mins"), 2), " minutes.\n", sep = "")
}


data.adults <- output[[1]] %>% 
  pivot_wider(., 
              id_cols=c(tox.scenario, simulation.ID, MTD.dose.level, sample.size), 
              names_from = dose.level, 
              values_from = c(doses, mean.tox, true.tox, toxicities.per.dose, patients.per.dose), 
              names_sep = "")

posterior <- output[[2]]


saveRDS(data.adults, "data_adults.rds")
data.adults <- readRDS("data_adults.rds")

saveRDS(posterior, "posterior.rds")
posterior <- readRDS("posterior.rds")

########################################
######## Set up pediatric trials
########################################


data <- tibble(data.adults)

# add pediatric doses to adult data 
data <- data %>% 
  mutate("MTD" = doses.adults[MTD.dose.level]) %>% 
  mutate("ped.dose1"=0.7*MTD, "ped.dose2"=MTD, "ped.dose3"=1.3*MTD, "ped.dose4"=1.6*MTD)


# add mean(s) of posterior to data
posterior <- tibble(posterior)
posterior.summary <- posterior %>% 
  group_by(tox.scenario,simulation.ID) %>%
  summarise(alpha = mean(V1), beta = mean(V2))
data <- left_join(data, posterior.summary, by=c("tox.scenario", "simulation.ID"))

# create skeleton
PLM2.mean.tox <- function(alpha, beta, di, dR){inv.logit(alpha + exp(beta) * log(di/dR))}
data <- data %>% 
  group_by(tox.scenario,simulation.ID) %>%
  do(data.frame("skeleton1"=PLM2.mean.tox(alpha=.$alpha, beta=.$beta, di=.$ped.dose1, dR=dR.adults),
                "skeleton2"=PLM2.mean.tox(alpha=.$alpha, beta=.$beta, di=.$ped.dose2, dR=dR.adults),
                "skeleton3"=PLM2.mean.tox(alpha=.$alpha, beta=.$beta, di=.$ped.dose3, dR=dR.adults),
                "skeleton4"=PLM2.mean.tox(alpha=.$alpha, beta=.$beta, di=.$ped.dose4, dR=dR.adults))) %>%
  left_join(data, ., by=c("tox.scenario", "simulation.ID"))

# full 1P-prior: add quantile of MTD according to posterior distribution and find the corresponding prior sd
quantile.MTD <- function(alpha,beta,MTD,dR){quantile(inv.logit(ifelse(alpha + exp(beta) * log(MTD/dR) > 8, 1, alpha + exp(beta) * log(MTD/dR))), probs=0.90)}
find.sd <- function(sd, skeleton2, target.quantile){abs(target.quantile-skeleton2^exp(qnorm(p=0.10, mean=0, sd=sd)))} 
data <- posterior %>% left_join(select(data, tox.scenario,simulation.ID,MTD), by=c("tox.scenario", "simulation.ID")) %>%   
  group_by(tox.scenario,simulation.ID) %>%
  do(data.frame("quantile.MTD"=quantile.MTD(alpha=.$V1, beta=.$V2, MTD=.$MTD, dR=dR.adults))) %>%
  left_join(data, ., by=c("tox.scenario", "simulation.ID")) 

specification.1P.full <- data %>% 
  group_by(tox.scenario,simulation.ID) %>%
  do(data.frame("prior.1P"=optimize(find.sd, skeleton2=.$skeleton2, target.quantile=.$quantile.MTD, interval=c(0.01, 2))$minimum)) 
specification.1P.full <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose")) %>% left_join(., specification.1P.full, by=c("tox.scenario", "simulation.ID")) %>% add_column(specification="1P.full")

# weak 1P-prior
specification.1P.weak <- data %>% 
  group_by(tox.scenario,simulation.ID) %>%
  do(data.frame("prior.1P"=optimize(find.sd, skeleton2=.$skeleton2, target.quantile=0.8, interval=c(0.01, 2))$minimum)) 
specification.1P.weak <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose")) %>% left_join(., specification.1P.weak, by=c("tox.scenario", "simulation.ID")) %>% add_column(specification="1P.weak")


########################################
######## Calibrate priors for 1P pediatric trials
########################################


##### prior calibration rules for first dose behavior
# if 2/2 tox, terminate trial
# if 1/2 tox, don't escalate
# if 2/4 tox, don't escalate
# if 0/2 tox, escalate

sd.calibrated <- list()

ncores <- detectCores()*3/4; cl <- makeCluster(ncores); registerDoParallel(cl);
start.time <- Sys.time()

for(t in 1:length(tox.adults)){
  
  specification.1P.full.t <- specification.1P.full %>% filter(tox.scenario==t)
  
  sd.calibrated <- foreach(s=1:n.sim, .combine=rbind, .multicombine=F, .packages = c("rstan", "dplyr"), .export = c("specification.1P.full.t", "crm.fixed")) %dopar% {
    
    condition <- rep(0, 7) # just starting value
    
    doses.current <- specification.1P.full.t %>% filter(simulation.ID==s & tox.scenario==t) %>% select(starts_with("ped.dose")) %>% as.matrix() %>% c()
    K <- length(doses.current)
    skeleton.current <- specification.1P.full.t %>% filter(simulation.ID==s & tox.scenario==t) %>% select(starts_with("skeleton")) %>% as.matrix() %>% c()
    mu <- 0
    sigma <- specification.1P.full.t %>% filter(simulation.ID==s & tox.scenario==t) %>% select(prior.1P) %>% as.matrix() %>% c()
    
    # there are seven conditions, and if condition is fulfilled, it is marked with '1' else '0'
    # hence the sum of all conditions should be 7
    
    while (sum(condition) < 7 & sigma<1.34){
      
      # if 2/2, recommended dose should be 1 with estimated probability of toxicity larger than 35%
      tox <- c(2,0,0,0)
      N <- c(2,0,0,0)
      data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton.current)
      samples <- sampling(crm.fixed, data = data.stan, warmup = 1000, iter = 3500, chains = 4, cores = 1, thin = 1, refresh = 0)
      samples <- extract(samples)
      est.tox <- skeleton.current^(exp(mean(samples$beta)))
      
      condition[1] <- ifelse(which.min(abs(est.tox - 0.25)) == 1, 1, 0) # recommended dose level should be 1
      condition[2] <- ifelse(est.tox[1] > 0.35, 1, 0) # estimated tox of first dose should be larger than 35%
      
      
      # if 1/2, recommended dose should be 1 with estimated probability of toxicity not larger than 35%
      tox <- c(1,0,0,0)
      N <- c(2,0,0,0)
      data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton.current)
      samples <- sampling(crm.fixed, data = data.stan, warmup = 1000, iter = 3500, chains = 4, cores = 1, thin = 1, refresh = 0)
      samples <- extract(samples)
      est.tox <- skeleton.current^(exp(mean(samples$beta)))
      
      condition[3] <- ifelse(which.min(abs(est.tox - 0.25)) == 1, 1, 0) # recommended dose level should be 1
      condition[4] <- ifelse(est.tox[1] > 0.35, 0, 1) # estimated tox of first dose should be not larger than 35%
      
      
      # if 2/4, recommended dose should be 1 with estimated probability of toxicity not larger than 35%
      # tox <- c(2,0,0,0)
      # N <- c(4,0,0,0)
      # data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton.current)
      # samples <- sampling(crm.fixed, data = data.stan, warmup = 1000, iter = 3500, chains = 4, cores = 1, thin = 1, refresh = 0)
      # samples <- extract(samples)
      # est.tox <- skeleton.current^(exp(mean(samples$beta)))
      # 
      # condition[5] <- ifelse(which.min(abs(est.tox - 0.25)) == 1, 1, 0) # recommended dose level should be 1
      # condition[6] <- ifelse(est.tox[1] > 0.35, 0, 1) # estimated tox of first dose should be not larger than 35%
      condition[5] <- 1
      condition[6] <- 1
      
      
      # if 0/2, recommended dose should be 1 with estimated probability of toxicity not larger than 35%
      tox <- c(0,0,0,0)
      N <- c(2,0,0,0)
      data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "mu" = mu, "sigma" = sigma, "skeleton" = skeleton.current)
      samples <- sampling(crm.fixed, data = data.stan, warmup = 1000, iter = 3500, chains = 4, cores = 1, thin = 1, refresh = 0)
      samples <- extract(samples)
      est.tox <- skeleton.current^(exp(mean(samples$beta)))
      
      condition[7] <- ifelse(which.min(abs(est.tox - 0.25)) >= 2, 1, 0) # recommended dose level should be 2 or larger
      
      sigma <- sigma + 0.05
    }
    return(c("sd.calibrated" = sigma-0.05))
  }
  
  sd.calibrated <- as.vector(sd.calibrated)
  if(t==1){data.prior.1P.calibrated <- tibble("tox.scenario"=t, "simulation.ID"=1:length(sd.calibrated), "prior.1P"=sd.calibrated)} else {
    data.prior.1P.calibrated <- rbind(data.prior.1P.calibrated, tibble("tox.scenario"=t, "simulation.ID"=1:length(sd.calibrated), "prior.1P"=sd.calibrated))
  }
}

stopCluster(cl); end.time <- Sys.time(); time.taken <- difftime(end.time, start.time, units = "min") # end parallel computing 
time.taken

saveRDS(data.prior.1P.calibrated, "data.prior.1P.calibrated.rds")
data.prior.1P.calibrated <- readRDS("data.prior.1P.calibrated.rds")

specification.1P.calibrated <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose")) %>% left_join(., data.prior.1P.calibrated, by=c("tox.scenario", "simulation.ID")) %>% add_column(specification="1P.calibrated")





############################################
########### prior for 2 parameter model
############################################



# add var.alpha, covariance and var.beta
specification.2P.full <- posterior %>% 
  group_by(tox.scenario,simulation.ID) %>%
  do(tibble("var.alpha" = c(cov(matrix(c(.$V1, .$V2), ncol=2, byrow=F)))[1],
            "covariance" = c(cov(matrix(c(.$V1, .$V2), ncol=2, byrow=F)))[2],
            "var.beta" = c(cov(matrix(c(.$V1, .$V2), ncol=2, byrow=F)))[4])) %>% 
  mutate("rho" = covariance/(sqrt(var.alpha) * sqrt(var.beta))) %>% 
  add_column(specification="2P.full", dR=dR.adults) %>% ungroup()
specification.2P.full <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose"), "alpha", "beta") %>% left_join(., specification.2P.full, by=c("tox.scenario", "simulation.ID"))


specification.2P.weak <- data %>% group_by(tox.scenario,simulation.ID) %>%
  do(tibble("alpha"=logit(0.25),
            "beta"=0.4606,
            "var.alpha"=1.939^2,
            "var.beta"=1.5^2,
            "covariance"=0,"dR"=.$ped.dose2)) %>% add_column(specification="2P.weak") %>% ungroup()
specification.2P.weak <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose")) %>% left_join(., specification.2P.weak, by=c("tox.scenario", "simulation.ID"))


#plm2.prior <- function(mu,var.alpha,var.beta,di,dR,target.quantile){abs(target.quantile - quantile(inv.logit(rnorm(1e4,mean=mu[1],sd=sqrt(var.alpha))+exp(rnorm(1e4,mean=mu[2],sd=sqrt(var.beta)))*log(di/dR)), probs=0.9))}
#data %>% group_by(tox.scenario,simulation.ID) %>% do(data.frame("var.alpha.weak"=optimize(plm2.prior, mu=c(.$alpha,.$beta), var.beta=1, di=.$ped.dose2, dR=dR.adults, target.quantile=0.8, interval=c(0.01, 9))$minimum)) %>% print(n=25)

# 
# # add the partial borrowing
# specification.2P.25 <- specification.2P.full %>% 
#   group_by(tox.scenario,simulation.ID) %>% 
#   do(data.frame("var.alpha"=.$var.alpha/0.25, "var.beta"=.$var.beta/0.25, "covariance"=.$rho * (sqrt(.$var.alpha/0.25) * sqrt(.$var.beta/0.25)))) %>% 
#   add_column(specification="2P.weighted.25", dR=dR.adults) %>% ungroup()
# specification.2P.25 <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose"), "alpha", "beta") %>% left_join(., specification.2P.25, by=c("tox.scenario", "simulation.ID"))
# 
# specification.2P.50 <- specification.2P.full %>% 
#   group_by(tox.scenario,simulation.ID) %>% 
#   do(data.frame("var.alpha"=.$var.alpha/0.50, "var.beta"=.$var.beta/0.50, "covariance"=.$rho * (sqrt(.$var.alpha/0.50) * sqrt(.$var.beta/0.50)))) %>% 
#   add_column(specification="2P.weighted.50", dR=dR.adults) %>% ungroup()
# specification.2P.50 <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose"), "alpha", "beta") %>% left_join(., specification.2P.50, by=c("tox.scenario", "simulation.ID"))
# 
# specification.2P.75 <- specification.2P.full %>% 
#   group_by(tox.scenario,simulation.ID) %>% 
#   do(data.frame("var.alpha"=.$var.alpha/0.75, "var.beta"=.$var.beta/0.75, "covariance"=.$rho * (sqrt(.$var.alpha/0.75) * sqrt(.$var.beta/0.75)))) %>% 
#   add_column(specification="2P.weighted.75", dR=dR.adults) %>% ungroup()
# specification.2P.75 <- data %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose"), "alpha", "beta") %>% left_join(., specification.2P.75, by=c("tox.scenario", "simulation.ID"))


########################################
######## Calibrate priors for 2P pediatric trials
########################################


ncores <- detectCores()*2/4; cl <- makeCluster(ncores); registerDoParallel(cl);
start.time <- Sys.time()

for(t in 1:length(tox.adults)){
  
  specification.2P.full.t <- specification.2P.full %>% filter(tox.scenario==t)
  
  
  Sigma.calibrated <- foreach(s=1:n.sim, .combine=rbind, .multicombine=F, .packages = c("rstan", "dplyr"), .export = c("specification.2P.full.t", "PLM2.fixed")) %dopar% {
    
    condition <- rep(0, 7) # just starting value
    
    doses.current <- specification.2P.full.t %>% filter(simulation.ID==s & tox.scenario==t) %>% select(starts_with("ped.dose")) %>% as.matrix() %>% c()
    K <- length(doses.current)
    data_Sigma <- specification.2P.full.t %>% filter(simulation.ID==s & tox.scenario==t)
    mu_alpha <- data_Sigma$alpha
    mu_beta <- data_Sigma$beta
    Sigma <- matrix(c(data_Sigma$var.alpha, data_Sigma$covariance, data_Sigma$covariance, data_Sigma$var.beta), ncol=2, byrow=T)
    rho <- data_Sigma$rho
    
    
    while (sum(condition) < 7 & Sigma[1,1]<prior.PLM2$Sigma[1,1]){
      
      # if 2/2, recommended dose should be 1 with estimated probability of toxicity larger than 53%
      tox <- c(2,0,0,0)
      N <- c(2,0,0,0)
      data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "dR"=dR.adults, "mu"=c(mu_alpha,mu_beta), "Sigma"=Sigma)
      samples <- sampling(PLM2.fixed, data = data.stan, warmup = 1000, iter = 11000, chains = 4, cores = 1, thin = 1, refresh = 0)
      samples <- rstan::extract(samples)
      est.tox <- inv.logit(mean(samples$logalpha) + exp(mean(samples$logbeta)) * log(doses.current/dR.adults))
      
      condition[1] <- ifelse(which.min(abs(est.tox - 0.25)) == 1, 1, 0) # recommended dose level should be 1
      condition[2] <- ifelse(est.tox[1] > 0.35, 1, 0) # estimated tox of first dose should be larger than 35%
      
      
      # if 1/2, recommended dose should be 1 with estimated probability of toxicity not larger than 35%
      tox <- c(1,0,0,0)
      N <- c(2,0,0,0)
      data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "dR"=dR.adults, "mu"=c(mu_alpha,mu_beta), "Sigma"=Sigma)
      samples <- sampling(PLM2.fixed, data = data.stan, warmup = 1000, iter = 11000, chains = 4, cores = 1, thin = 1, refresh = 0)
      samples <- rstan::extract(samples)
      est.tox <- inv.logit(mean(samples$logalpha) + exp(mean(samples$logbeta)) * log(doses.current/dR.adults))
      
      condition[3] <- ifelse(which.min(abs(est.tox - 0.25)) == 1, 1, 0) # recommended dose level should be 1
      condition[4] <- ifelse(est.tox[1] > 0.35, 0, 1) # estimated tox of first dose should be not larger than 35%
      
      
      # if 2/4, recommended dose should be 1 with estimated probability of toxicity not larger than 50%
      # tox <- c(2,0,0,0)
      # N <- c(4,0,0,0)
      # data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "dR"=dR.adults, "mu"=c(mu_alpha,mu_beta), "Sigma"=Sigma)
      # samples <- sampling(PLM2.fixed, data = data.stan, warmup = 1000, iter = 11000, chains = 4, cores = 1, thin = 1, refresh = 0)
      # samples <- rstan::extract(samples)
      # est.tox <- inv.logit(mean(samples$logalpha) + exp(mean(samples$logbeta)) * log(doses.current/dR.adults))
      # 
      # condition[5] <- ifelse(which.min(abs(est.tox - 0.25)) == 1, 1, 0) # recommended dose level should be 1
      # condition[6] <- ifelse(est.tox[1] > 0.35, 0, 1) # estimated tox of first dose should be not larger than 35%
      condition[5] <- 1
      condition[6] <- 1
      
      # if 0/2, recommendation should be dose >= 2
      tox <- c(0,0,0,0)
      N <- c(2,0,0,0)
      data.stan <- list("K"=K, "tox"=tox, "N"=N, "doses"=doses.current, "dR"=dR.adults, "mu"=c(mu_alpha,mu_beta), "Sigma"=Sigma)
      samples <- sampling(PLM2.fixed, data = data.stan, warmup = 1000, iter = 11000, chains = 4, cores = 1, thin = 1, refresh = 0)
      samples <- rstan::extract(samples)
      est.tox <- inv.logit(mean(samples$logalpha) + exp(mean(samples$logbeta)) * log(doses.current/dR.adults))
      
      condition[7] <- ifelse(which.min(abs(est.tox - 0.25)) >= 2, 1, 0) # recommended dose level should 2 or larger
      
      
      Sigma[1,1] <- Sigma[1,1] + 0.25
      Sigma[2,2] <- Sigma[2,2] + 0.1
      Sigma[2,1] <- sqrt(Sigma[1,1])*sqrt(Sigma[2,2])*rho
      Sigma[1,2] <- sqrt(Sigma[1,1])*sqrt(Sigma[2,2])*rho
    }
    Sigma[1,1] <- Sigma[1,1] - 0.25
    Sigma[2,2] <- Sigma[2,2] - 0.1 
    Sigma[2,1] <- sqrt(Sigma[1,1])*sqrt(Sigma[2,2])*rho
    Sigma[1,2] <- sqrt(Sigma[1,1])*sqrt(Sigma[2,2])*rho
    
    return(c("Sigma.calibrated" = Sigma, "rho"=rho))
  }
  Sigma.calibrated <- as.matrix(Sigma.calibrated)
  if(t==1){data.prior.2P.calibrated <- tibble("tox.scenario"=t, "simulation.ID"=1:nrow(Sigma.calibrated), "var.alpha"=Sigma.calibrated[,1], "var.beta"=Sigma.calibrated[,4], "covariance"=Sigma.calibrated[,2], "rho"=Sigma.calibrated[,5], "dR"=dR.adults)
  } else {
    data.prior.2P.calibrated <- rbind(data.prior.2P.calibrated, tibble("tox.scenario"=t, "simulation.ID"=1:nrow(Sigma.calibrated), "var.alpha"=Sigma.calibrated[,1], "var.beta"=Sigma.calibrated[,4], "covariance"=Sigma.calibrated[,2], "rho"=Sigma.calibrated[,5], "dR"=dR.adults))
  }
}

stopCluster(cl); end.time <- Sys.time(); time.taken <- difftime(end.time, start.time, units = "min") # end parallel computing 
time.taken

saveRDS(data.prior.2P.calibrated, "data.prior.2P.calibrated.rds")
data.prior.2P.calibrated <- readRDS("data.prior.2P.calibrated.rds")
specification.2P.calibrated <- specification.2P.full %>% select("tox.scenario", "simulation.ID", starts_with("skeleton"), starts_with("ped.dose"), "alpha", "beta") %>% left_join(., data.prior.2P.calibrated, by=c("tox.scenario", "simulation.ID")) %>% add_column(specification="2P.calibrated")



#######################################
# mixture priors
specification.1P.mixture <- cbind(subset(specification.1P.full, select=-c(specification)), "prior.1P2"=specification.1P.weak$prior.1P) %>% add_column(specification="1P.mixture")

specification.2P.mixture <- cbind(subset(specification.2P.full, select=-c(specification)), "var.alpha2"=prior.PLM2$Sigma[1,1], "var.beta2"=prior.PLM2$Sigma[2,2], "covariance2"=0) %>% add_column(specification="2P.mixture")

# combine all prior inputs together
full.specifications <- bind_rows(specification.1P.weak,specification.1P.full,specification.1P.mixture,specification.1P.calibrated,
                                 specification.2P.weak,specification.2P.full,specification.2P.mixture,specification.2P.calibrated,
                                 #specification.2P.25,specification.2P.50,specification.2P.75
) %>% 
  add_column(sample.size=10, ped.tox.weaker="weaker", ped.tox.same="same", ped.tox.stronger="stronger") %>% 
  pivot_longer(., cols=starts_with("ped.tox."), values_to="ped.scenario") %>% select(-c(name))

full.specifications <- full.specifications %>% 
  rowwise() %>% 
  mutate(tox.dose1 = case_when(ped.scenario=="weaker" ~ fun.true.tox.ped.weaker(doses=ped.dose1,t=tox.scenario),
                               ped.scenario=="same" ~ fun.true.tox.ped.same(doses=ped.dose1,t=tox.scenario),
                               ped.scenario=="stronger" ~ fun.true.tox.ped.stronger(doses=ped.dose1,t=tox.scenario)),
         tox.dose2 = case_when(ped.scenario=="weaker" ~ fun.true.tox.ped.weaker(doses=ped.dose2,t=tox.scenario),
                               ped.scenario=="same" ~ fun.true.tox.ped.same(doses=ped.dose2,t=tox.scenario),
                               ped.scenario=="stronger" ~ fun.true.tox.ped.stronger(doses=ped.dose2,t=tox.scenario)),
         tox.dose3 = case_when(ped.scenario=="weaker" ~ fun.true.tox.ped.weaker(doses=ped.dose3,t=tox.scenario),
                               ped.scenario=="same" ~ fun.true.tox.ped.same(doses=ped.dose3,t=tox.scenario),
                               ped.scenario=="stronger" ~ fun.true.tox.ped.stronger(doses=ped.dose3,t=tox.scenario)),
         tox.dose4 = case_when(ped.scenario=="weaker" ~ fun.true.tox.ped.weaker(doses=ped.dose4,t=tox.scenario),
                               ped.scenario=="same" ~ fun.true.tox.ped.same(doses=ped.dose4,t=tox.scenario),
                               ped.scenario=="stronger" ~ fun.true.tox.ped.stronger(doses=ped.dose4,t=tox.scenario))) %>% ungroup()



################################################
# oracle models

# 1P-CRM: skeleton matches exactly the toxicity probabilities with same variance as under full borrowing
full.specifications <- rbind(full.specifications,
                             full.specifications %>% filter(specification=="1P.full") %>% 
                               mutate(skeleton1=tox.dose1,skeleton2=tox.dose2,skeleton3=tox.dose3,skeleton4=tox.dose4, specification="1P.oracle"))


# 2P-CRM: dR is the dose that actually has 30% tox prob, di is the dose with 10% tox prob, Sigma rather uninformative (2,0,0,1)

# the following is the same for all simulations within one scenario:
# find doses with 10% and 30% tox prob and make them dR and di
# find corresponding alpha and beta
# make the variance of oracle approach be the same as full borrowing

# df.comb.scenarios <- tibble(expand.grid("tox.scenario"=unique(full.specifications$tox.scenario), "ped.scenario"=unique(full.specifications$ped.scenario))) %>%
#   rowwise() %>%
#   mutate(dose.di = case_when(ped.scenario=="weaker" ~ optimize(function(doses){abs(fun.true.tox.ped.weaker(doses=doses,t=tox.scenario)-0.1)},1:300)$minimum,
#                              ped.scenario=="same" ~ optimize(function(doses){abs(fun.true.tox.ped.same(doses=doses,t=tox.scenario)-0.1)},1:300)$minimum,
#                              ped.scenario=="stronger" ~ optimize(function(doses){abs(fun.true.tox.ped.stronger(doses=doses,t=tox.scenario)-0.1)},1:300)$minimum),
#          dR = case_when(ped.scenario=="weaker" ~ round(optimize(function(doses){abs(fun.true.tox.ped.weaker(doses=doses,t=tox.scenario)-0.3)},1:300)$minimum),
#                              ped.scenario=="same" ~ round(optimize(function(doses){abs(fun.true.tox.ped.same(doses=doses,t=tox.scenario)-0.3)},1:300)$minimum),
#                              ped.scenario=="stronger" ~ round(optimize(function(doses){abs(fun.true.tox.ped.stronger(doses=doses,t=tox.scenario)-0.3)},1:300)$minimum)),
#          alpha = logit(0.3), beta = (logit(0.1)-logit(0.3))/log(dose.di/dR)) %>% ungroup() %>% select(-c(dose.di))
# 
# full.specifications <- bind_rows(full.specifications, 
#                                  full.specifications %>% filter(specification=="2P.full") %>% select(-c(dR,alpha,beta,specification)) %>% 
#                                    left_join(.,df.comb.scenarios,by=c("tox.scenario", "ped.scenario")) %>% add_column(specification="2P.oracle"))

oracle.2P <- function(){
  dose.data <- tibble("dose"=5:80, "tox"=tox.fun(5:80)) %>% filter(tox>0, tox<1)
  dose.data2 <- dose.data %>% filter(tox>0.15, tox<0.35)
  # we want the intercept to be the MTD
  dR <- dose.data2$dose[which(abs(dose.data2$tox-0.25) == min(abs(dose.data2$tox-0.25)))]
  dose.data2$x <- dose.data2$dose/dR
  mod <- lm(I(qlogis(tox)-qlogis(tox.fun(dR))) ~ 0 + I(log(dose/dR)), data=dose.data2)
  return(data.frame("alpha"=qlogis(tox.fun(dR)), "beta"=coef(mod), "dR"=dR))
}
tox.fun <- function(doses){dose.tox.weak(doses)}
oracle.2P()


tox.fun <- function(doses){dose.tox.weak(doses)};
d1 <- tibble(tox.scenario=1,ped.scenario="same",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){dose.tox.weak(doses)^1.3}; 
d2 <- tibble(tox.scenario=1,ped.scenario="weaker",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){dose.tox.weak(doses)^0.7}
d3 <- tibble(tox.scenario=1,ped.scenario="stronger",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)

tox.fun <- function(doses){dose.tox.moderate(doses)}
d4 <- tibble(tox.scenario=2,ped.scenario="same",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){dose.tox.moderate(doses)^1.3}
d5 <- tibble(tox.scenario=2,ped.scenario="weaker",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){dose.tox.moderate(doses)^0.7}
d6 <- tibble(tox.scenario=2,ped.scenario="stronger",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)

tox.fun <- function(doses){dose.tox.strong(doses)}
d7 <- tibble(tox.scenario=3,ped.scenario="same",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){dose.tox.strong(doses)^1.3}
d8 <- tibble(tox.scenario=3,ped.scenario="weaker",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){dose.tox.strong(doses)^0.4}
d9 <- tibble(tox.scenario=3,ped.scenario="stronger",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)

tox.fun <- function(doses){spline.plateau(doses)}
d10 <- tibble(tox.scenario=4,ped.scenario="same",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){spline.plateau.weaker(doses)}
d11 <- tibble(tox.scenario=4,ped.scenario="weaker",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){spline.plateau.stronger(doses)}
d12 <- tibble(tox.scenario=4,ped.scenario="stronger",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)

par(mfrow=c(1,3))
tox.fun <- function(doses){spline.waves(doses)}
d13 <- tibble(tox.scenario=5,ped.scenario="same",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){spline.waves.weaker(doses)}
d14 <- tibble(tox.scenario=5,ped.scenario="weaker",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)
tox.fun <- function(doses){spline.waves.stronger(doses)}
d15 <- tibble(tox.scenario=5,ped.scenario="stronger",dR=oracle.2P()$dR,alpha=oracle.2P()$alpha,beta=oracle.2P()$beta)

data.oracle.2P <- cbind(rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15),"specification"="2P.oracle")

full.specifications <-  rbind(full.specifications,
                              full.specifications %>% filter(specification=="2P.full") %>% select(-c(dR,alpha,beta,specification)) %>% 
                                left_join(.,data.oracle.2P,by=c("tox.scenario", "ped.scenario")))


#for quicker simulation:
#full.specifications <- full.specifications %>% filter(simulation.ID<=100)
# full.specifications <- full.specifications %>% filter(specification=="2P.oracle", tox.scenario %in% c(4,5))


full.specifications$ped.id <- as.numeric(rownames(full.specifications))

saveRDS(full.specifications, "full_specifications.rds")
full.specifications <- readRDS("full_specifications.rds")


###### run parallelized simulation
# storing of simulation data
results <- list()


# start parallelization
ncores <- detectCores()*3/4; cl <- makeCluster(ncores); registerDoParallel(cl);

steps <- 250 # in how many steps will the simulation be partitioned? One step should include not more than approx. 1000 trials to avoid running out of memory.

for(i in 1:steps){
  
  start.time <- Sys.time()
  
  results <- foreach(index=((nrow(full.specifications)/steps)*(i-1)+1):((nrow(full.specifications)/steps)*i), .combine=rbind, .multicombine=F, .packages = c("rstan"),
                     .export = c("full.specifications","sim.phaseI","crm.fixed","PLM2.fixed","crm.mixture","PLM2.mixture")) %dopar% {
                       
                       # define columns that are used for the prior
                       input.simulation <- full.specifications[index,]
                       
                       if(input.simulation$specification %in% c("1P.weak", "1P.full", "1P.calibrated", "1P.oracle")){
                         prior <- list(0,input.simulation$prior.1P)
                         alg <- "CRM"; borrow <- "fixed"; stanmodel <- crm.fixed
                       } else if(input.simulation$specification=="1P.mixture"){
                         prior <- list(c(0,input.simulation$prior.1P), c(0,input.simulation$prior.1P2))
                         alg <- "CRM"; borrow <- "mixture"; stanmodel <- crm.mixture
                       } else if(input.simulation$specification %in% c("2P.weak", "2P.full", "2P.weighted", "2P.calibrated", "2P.weighted.25", "2P.weighted.50", "2P.weighted.75", "2P.oracle")){
                         prior <- list(input.simulation$alpha,input.simulation$beta, matrix(c(input.simulation$var.alpha,rep(input.simulation$covariance,2),input.simulation$var.beta),nrow=2,byrow=T))
                         alg <- "PLM2"; borrow <- "fixed"; stanmodel <- PLM2.fixed
                       } else if(input.simulation$specification=="2P.mixture"){
                         prior <- list(list(c(input.simulation$alpha,input.simulation$beta), matrix(c(input.simulation$var.alpha,rep(input.simulation$covariance,2),input.simulation$var.beta),nrow=2,byrow=T)),
                                       list(c(input.simulation$alpha,input.simulation$beta), matrix(c(input.simulation$var.alpha2,rep(input.simulation$covariance2,2),input.simulation$var.beta2),nrow=2,byrow=T)))
                         alg <- "PLM2"; borrow <- "mixture"; stanmodel <- PLM2.mixture
                       }
                       
                       results.sim <- sim.phaseI(algorithm=alg, doses = c(input.simulation$ped.dose1,input.simulation$ped.dose2,input.simulation$ped.dose3,input.simulation$ped.dose4), 
                                                 target.tox = 0.25, true.tox = c(input.simulation$tox.dose1,input.simulation$tox.dose2,input.simulation$tox.dose3,input.simulation$tox.dose4), 
                                                 prior = prior, borrow = borrow, prior.distr = "Normal",
                                                 skeleton = c(input.simulation$skeleton1,input.simulation$skeleton2,input.simulation$skeleton3,input.simulation$skeleton4), 
                                                 dR = input.simulation$dR, stopping.rule = "mean", escalation.rule = "mean",
                                                 sample.size = input.simulation$sample.size, n.sim = 1, cohort.size = 2, 
                                                 iterations = 5000, stanmodel = stanmodel)
                       
                       results.sim <- cbind(results.sim[[1]], "index"=index)                                              
                       return(results.sim)
                     }
  
  if(i==1) { ped.data <- results } else { ped.data <- rbind(ped.data, results) }
  
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units = "min") 
  print(paste("Iteration Number", (nrow(full.specifications)/steps)*i, ", time taken for this iteration", time.taken, ", current time", Sys.time()))
}

stopCluster(cl) # end parallel computing 

saveRDS(ped.data, "ped_data.rds")
ped.data <- readRDS("ped_data.rds")



