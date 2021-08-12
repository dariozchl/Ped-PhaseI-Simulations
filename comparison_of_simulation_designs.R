###############
### compare simulation designs with various number of pediatric trials subsequent to a adult trial: 
### 250*40 vs 1e4 * 1
###############

library(rstan)
library(doParallel)
library(tidyverse)

PLM2.fixed <- stan_model("BLRM.stan")
PLM2.mixture <- stan_model("BLRM_mixture.stan")

source("crm_function.R")

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


################### toxicity scenarios ###################

LL4 <- function(x,b,c,d,e){c + (d-c)/(1+exp(b*(log(x)-log(e))))}
x=seq(1,100,by=0.1)
dose.tox.moderate <- function(doses){1-LL4(x=doses,b=1.3,c=0,d=1,e=50)}
tox.adults.moderate <- dose.tox.moderate(doses.adults)

ped.tox.stronger <- function(doses){dose.tox.moderate(doses)^0.7}
ped.tox.same <- function(doses){dose.tox.moderate(doses)}


################### Adult Prior ###################
# reference dose
dR.adults <- 20

# specify weakly informative prior for 2P-CRM
prior.PLM2 <- list(mu=c(logit(0.25), 0.4606), Sigma=matrix(c(1.939^2, 0, 0, 1.5^2), nrow = 2, byrow = T))



########################################
######## Simulation adult trials
########################################

start.time <- Sys.time()

n.sim <- 5e4

library(doParallel)
library(tidyverse)
cores <- detectCores()*3/4
cl <- makeCluster(cores)
registerDoParallel(cl)

sim_results <- foreach(nsim = 1:n.sim, .packages=c("rstan", "tidyverse")) %dopar% { 
  index <- nsim
  output.temporary <- NULL
  while(!any(output.temporary[[1]]$mean.tox > (1/6) & output.temporary[[1]]$mean.tox < (1/3))){
    output.temporary <- sim.phaseI(algorithm="PLM2", doses = doses.adults, target.tox = 0.25, 
                                   true.tox = tox.adults.moderate,
                                   prior = list(unlist(prior.PLM2$mu)[1], unlist(prior.PLM2$mu)[2], prior.PLM2$Sigma), 
                                   borrow = "fixed", prior.distr = "Normal",
                                   skeleton = NULL, dR = dR.adults, sample.size = ss.adults, cohort.size = 3,
                                   n.sim = 1, iterations = 10000, stanmodel = PLM2.fixed,
                                   stopping.rule = "mean", escalation.rule = "mean")
  } 
  
  return(list(
    #simulation results:
    as_tibble(output.temporary[[1]]) %>% 
      pivot_wider(., 
                  id_cols=c(simulation.ID, MTD.dose.level, sample.size), 
                  names_from = dose.level, 
                  values_from = c(doses, mean.tox, true.tox, toxicities.per.dose, patients.per.dose), 
                  names_sep = "") %>% mutate("simulation.ID"=nsim), 
    #posterior:
    as_tibble(output.temporary[[2]], .name_repair = "universal") %>% 
      do(tibble("alpha" = mean(.$...1), "beta" = mean(.$...2),
                "var.alpha" = c(cov(matrix(c(.$...1, .$...2), ncol=2, byrow=F)))[1],
                "covariance" = c(cov(matrix(c(.$...1, .$...2), ncol=2, byrow=F)))[2],
                "var.beta" = c(cov(matrix(c(.$...1, .$...2), ncol=2, byrow=F)))[4])) %>% 
      mutate("rho" = covariance/(sqrt(var.alpha) * sqrt(var.beta)),
             "simulation.ID"=nsim))
  )
}

stopCluster(cl)
Sys.time() - start.time

saveRDS(sim_results, "comparison_of_simulation_designs_adult_data.RDS")
sim_results <- readRDS("comparison_of_simulation_designs_adult_data.RDS")


# data preparation
sim_results_tibbles <- lapply(purrr::transpose(sim_results), function(lst) do.call(rbind, lst))
data <- sim_results_tibbles[[1]]
posterior.summary <- sim_results_tibbles[[2]]


# add pediatric doses to adult data 
data <- data %>% 
  mutate("MTD" = doses.adults[MTD.dose.level]) %>% 
  mutate("ped.dose1"=0.7*MTD, "ped.dose2"=MTD, "ped.dose3"=1.3*MTD, "ped.dose4"=1.6*MTD)


########################################
######## Prior derivation
########################################

# add mean(s) of posterior to data
data <- left_join(data, posterior.summary, by=c("simulation.ID"))

specification.2P.full <- data %>% select("simulation.ID", starts_with("ped.dose"), "alpha", "beta", "covariance", starts_with("var"), "rho") %>% add_column(specification="2P.full", "dR"=dR.adults)


specification.2P.weak <- data %>% select("simulation.ID", starts_with("ped.dose")) %>%
  mutate("alpha"=logit(0.25),
         "beta"=0.4606,
         "var.alpha"=1.939^2,
         "var.beta"=1.5^2,
         "covariance"=0,"dR"=.$ped.dose2) %>% 
  add_column(specification="2P.weak")


specification.2P.mixture <- as_tibble(cbind(subset(specification.2P.full, select=-c(specification)), "var.alpha2"=prior.PLM2$Sigma[1,1], "var.beta2"=prior.PLM2$Sigma[2,2], "covariance2"=0) %>% 
                                        add_column(specification="2P.mixture"))


full.specifications <- bind_rows(specification.2P.weak %>% add_column(ped.scenario="same"), 
                                 specification.2P.mixture %>% add_column(ped.scenario="same"),
                                 specification.2P.weak %>% add_column(ped.scenario="stronger"), 
                                 specification.2P.mixture %>% add_column(ped.scenario="stronger")) 

full.specifications <- full.specifications %>% 
  rowwise() %>% 
  mutate(tox.dose1 = case_when(ped.scenario=="stronger" ~ ped.tox.stronger(doses=ped.dose1),
                               ped.scenario=="same" ~ ped.tox.same(doses=ped.dose1)),
         tox.dose2 = case_when(ped.scenario=="stronger" ~ ped.tox.stronger(doses=ped.dose2),
                               ped.scenario=="same" ~ ped.tox.same(doses=ped.dose2)),
         tox.dose3 = case_when(ped.scenario=="stronger" ~ ped.tox.stronger(doses=ped.dose3),
                               ped.scenario=="same" ~ ped.tox.same(doses=ped.dose3)),
         tox.dose4 = case_when(ped.scenario=="stronger" ~ ped.tox.stronger(doses=ped.dose4),
                               ped.scenario=="same" ~ ped.tox.same(doses=ped.dose4)),
         sample.size = 10)

########################################
######## Simulation pediatric trials 1x1
########################################

full.specifications$ped.id <- as.numeric(rownames(full.specifications))

ncores <- detectCores()*3/4; cl <- makeCluster(ncores); registerDoParallel(cl);

steps <- n.sim/1000 # in how many steps will the simulation be partitioned? One step should include not more than approx. 1000 trials to avoid running out of memory.

for(i in 1:steps){
  
  start.time <- Sys.time()
  
  results <- foreach(index=((nrow(full.specifications)/steps)*(i-1)+1):((nrow(full.specifications)/steps)*i), .combine=rbind, .multicombine=F, .packages = c("rstan")) %dopar% {
    
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
  
  if(i==1) { ped_data1x1 <- results } else { ped_data1x1 <- rbind(ped_data1x1, results) }
  
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units = "min") 
  print(paste("Iteration Number", (nrow(full.specifications)/steps)*i, ", time taken for this iteration", time.taken, ", current time", Sys.time()))
}

stopCluster(cl) # end parallel computing 

saveRDS(ped_data1x1, "ped_data1x1.rds")
ped_data1x1 <- readRDS("ped_data1x1.rds")



########################################
######## Simulation pediatric trials manyx1
########################################

sample.ID <- sample(unique(full.specifications$simulation.ID),n.sim/100,replace=FALSE)
full.specifications.manyx1 <- full.specifications %>% filter(simulation.ID %in% sample.ID) %>%
  rownames_to_column(var="index") %>% mutate(index = as.numeric(index))

ncores <- detectCores()*3/4; cl <- makeCluster(ncores); registerDoParallel(cl);

steps <- n.sim/1000 # in how many steps will the simulation be partitioned? One step should include not more than approx. 1000 trials to avoid running out of memory.

for(i in 1:steps){
  
  start.time <- Sys.time()
  
  results <- foreach(index=((nrow(full.specifications.manyx1)/steps)*(i-1)+1):((nrow(full.specifications.manyx1)/steps)*i), .combine=rbind, .multicombine=F, .packages = c("rstan")) %dopar% {
    
    # define columns that are used for the prior
    input.simulation <- full.specifications.manyx1[index,]
    
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
                              sample.size = input.simulation$sample.size, n.sim = 100, cohort.size = 2, 
                              iterations = 5000, stanmodel = stanmodel)
    
    results.sim <- cbind(results.sim[[1]], "index"=index)                                              
    return(results.sim)
  }
  
  if(i==1) { ped_datamanyx1 <- results } else { ped_datamanyx1 <- rbind(ped_datamanyx1, results) }
  
  end.time <- Sys.time()
  time.taken <- difftime(end.time, start.time, units = "min") 
  print(paste("Iteration Number", (nrow(full.specifications.manyx1)/steps)*i, ", time taken for this iteration", time.taken, ", current time", Sys.time()))
}

stopCluster(cl) # end parallel computing 

saveRDS(ped_datamanyx1, "ped_datamanyx1.rds")



########################################
######## Results pediatric trials 1 x 1
########################################

ped_data1x1 <- readRDS("ped_data1x1.rds")

ped_data1x1 <- as_tibble(ped_data1x1) %>% 
  mutate(toxicity.dose = case_when(true.tox>(1/3) ~ "overdosing",true.tox<(1/6) ~ "underdosing",TRUE ~ "acceptable")) %>% 
  group_by(simulation.ID, index) %>% 
  mutate(toxicity.MTD = case_when(
    MTD.dose.level==0 ~ "early ET",
    MTD.dose.level==1 ~ toxicity.dose[which(dose.level==1)], 
    MTD.dose.level==2 ~ toxicity.dose[which(dose.level==2)],
    MTD.dose.level==3 ~ toxicity.dose[which(dose.level==3)], 
    MTD.dose.level==4 ~ toxicity.dose[which(dose.level==4)])
  ) %>% 
  pivot_wider(., id_cols=c(simulation.ID, index, sample.size, MTD.dose.level, toxicity.MTD), names_from = dose.level, 
              values_from = c(doses, mean.tox, true.tox, toxicities.per.dose, patients.per.dose, toxicity.dose), 
              names_sep = "") %>% 
  mutate(toxicity.MTD = case_when(
    MTD.dose.level==0 & true.tox1>(1/3) ~ "correct ET", 
    MTD.dose.level==4 & true.tox4<(1/6) & mean.tox4<(1/6) ~ "correct ET", 
    MTD.dose.level==0 & true.tox1<(1/3)  ~ "false ET",
    TRUE ~ toxicity.MTD))

ped_data1x1 <- ped_data1x1 %>% mutate(ped.id = index)

ped_data1x1 <- left_join(ped_data1x1, full.specifications, by=c("ped.id"))


########################################
######## Results pediatric trials many x 1
########################################

ped_datamanyx1 <- readRDS("ped_datamanyx1.rds")

ped_datamanyx1 <- as_tibble(ped_datamanyx1) %>% 
  mutate(toxicity.dose = case_when(true.tox>(1/3) ~ "overdosing",true.tox<(1/6) ~ "underdosing",TRUE ~ "acceptable")) %>% 
  group_by(simulation.ID, index) %>% 
  mutate(toxicity.MTD = case_when(
    MTD.dose.level==0 ~ "early ET",
    MTD.dose.level==1 ~ toxicity.dose[which(dose.level==1)], 
    MTD.dose.level==2 ~ toxicity.dose[which(dose.level==2)],
    MTD.dose.level==3 ~ toxicity.dose[which(dose.level==3)], 
    MTD.dose.level==4 ~ toxicity.dose[which(dose.level==4)])
  ) %>% 
  pivot_wider(., id_cols=c(simulation.ID, index, sample.size, MTD.dose.level, toxicity.MTD), names_from = dose.level, 
              values_from = c(doses, mean.tox, true.tox, toxicities.per.dose, patients.per.dose, toxicity.dose), 
              names_sep = "") %>% 
  mutate(toxicity.MTD = case_when(
    MTD.dose.level==0 & true.tox1>(1/3) ~ "correct ET", 
    MTD.dose.level==4 & true.tox4<(1/6) & mean.tox4<(1/6) ~ "correct ET", 
    MTD.dose.level==0 & true.tox1<(1/3)  ~ "false ET",
    TRUE ~ toxicity.MTD))

ped_datamanyx1 <- ped_datamanyx1 %>% mutate(ped.id = full.specifications.manyx1[index,]$ped.id) %>% rename(simulation.ID.ped = simulation.ID)

ped_datamanyx1 <- left_join(ped_datamanyx1, full.specifications.manyx1, by=c("index"))



bind_rows(ped_datamanyx1 %>% add_column(design="many_by_1"), ped_data1x1 %>% add_column(design="1_by_1")) %>% 
  group_by("Prior" = unlist(specification), ped.scenario, design) %>% 
  summarize("Acceptable" = sum(toxicity.MTD=="acceptable"),
            "Correct ET" = sum(toxicity.MTD=="correct ET"),
            "Sum of \ncorrect decisions" = sum(toxicity.MTD=="acceptable" | toxicity.MTD=="correct ET"),
            "N"=n()) %>% 
  pivot_longer(., cols=c(Acceptable, `Correct ET`, `Sum of \ncorrect decisions`), names_to="Decision", values_to="Number") %>%
  rowwise() %>%
  mutate("Proportion" = 100*Number/N) %>%
  pivot_wider(id_cols = c(Prior, ped.scenario, design, Decision), 
              values_from=c(Proportion), 
              names_from=design, names_glue = '{design}.{.value}') %>%
  mutate(prop.diff = `1_by_1.Proportion` - many_by_1.Proportion)
