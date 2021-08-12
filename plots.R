
# Author: Dario Zocholl, Institute of Biometry and Clinical Epidemiology, Charité - Universitätsmedizin Berlin

########################################
######## This script performs comparison between simulation designs: 
#           - 1 x 1 (i.e. each adult trial is followed by one simulated pediatric trial) versus
#           - 1 x many (i.e. each adult trial is followed by many simulated pediatric trials)
######## All files are publicly available under https://github.com/dariozchl/Ped-PhaseI-Simulations
########################################


library(tidyverse)
library(ggpubr)
library(viridis)



inv.logit <- function(x){exp(x)/(1+exp(x))}
logit <- function(p){log(p/(1-p))}
ss.adults <- 30
doses.adults <- c(5,10,15,20,27,36,47)

################################
# toxicity scenarios

LL4 <- function(x,b,c,d,e){c + (d-c)/(1+exp(b*(log(x)-log(e))))}
x=seq(1,100,by=0.1)
dose.tox.weak <- function(doses){1-LL4(x=doses,b=1.5,c=0,d=1,e=70)}
dose.tox.moderate <- function(doses){1-LL4(x=doses,b=1.3,c=0,d=1,e=50)}
dose.tox.strong <- function(doses){1-LL4(x=doses,b=1.8,c=0,d=1,e=27)}
probs1 <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)
spline.plateau <- splinefun(c(3.5, doses.adults,188,40), probs1, "monoH.FC")
probs2 <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)^0.7
spline.plateau.stronger <- splinefun(c(3.5, doses.adults,188,40), probs2, "monoH.FC")
probs3 <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)^1.5
spline.plateau.weaker <- splinefun(c(3.5, doses.adults,188,40), probs3, "monoH.FC")
probs1 <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)
spline.waves <- splinefun(c(doses.adults,188), probs1, "monoH.FC")
probs2 <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)^0.7
spline.waves.stronger <- splinefun(c(doses.adults,188), probs2, "monoH.FC")
probs3 <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)^1.5
spline.waves.weaker <- splinefun(c(doses.adults,188), probs3, "monoH.FC")
tox.adults.weak <- dose.tox.weak(doses.adults)
tox.adults.moderate <- dose.tox.moderate(doses.adults)
tox.adults.strong <- dose.tox.strong(doses.adults)
tox.adults.spline.plateau <- spline.plateau(doses.adults)
tox.adults.spline.waves <- spline.waves(doses.adults)
tox.adults <- list(tox.adults.weak, tox.adults.moderate, tox.adults.strong, tox.adults.spline.plateau, tox.adults.spline.waves)


ped.data <- tibble(readRDS("ped_data.rds"))
full.specifications <- tibble(readRDS("full_specifications.rds"))

ped.data <- ped.data %>% mutate(toxicity.dose = case_when(true.tox>(1/3) ~ "overdosing",true.tox<(1/6) ~ "underdosing",TRUE ~ "acceptable"))

ped.data <- ped.data %>% group_by(index) %>% mutate(toxicity.MTD = case_when(MTD.dose.level==0 ~ "early ET",
                                                                             MTD.dose.level==1 ~ toxicity.dose[which(dose.level==1)], MTD.dose.level==2 ~ toxicity.dose[which(dose.level==2)],
                                                                             MTD.dose.level==3 ~ toxicity.dose[which(dose.level==3)], MTD.dose.level==4 ~ toxicity.dose[which(dose.level==4)]))

data <- ped.data %>% select(-simulation.ID) %>% pivot_wider(., id_cols=c(index, sample.size, MTD.dose.level, toxicity.MTD), names_from = dose.level, values_from = c(doses, mean.tox, true.tox, toxicities.per.dose, patients.per.dose, toxicity.dose), names_sep = "")

data <- data %>% mutate(toxicity.MTD = case_when(MTD.dose.level==0 & true.tox1>(1/3) ~ "correct ET", 
                                                 MTD.dose.level==4 & true.tox4<(1/6) & mean.tox4<(1/6) ~ "correct ET", 
                                                 MTD.dose.level==0 & true.tox1<(1/3)  ~ "false ET",
                                                 TRUE ~ toxicity.MTD)) %>% rename(ped.id = index)

data <- left_join(data, full.specifications, by="ped.id")

data <- data %>% 
  mutate(adult.scenario = case_when(tox.scenario==1 ~ "weak",tox.scenario==2 ~ "moderate",tox.scenario==3 ~ "strong",tox.scenario==4 ~ "spline.waves",tox.scenario==5 ~ "spline.plateau")) %>%
  rename(sample.size = sample.size.x)







################################################################################################################################################################################################################################################################################################################################################
# scenarios

doses.ped <- c(doses.adults*0.7,doses.adults*1,doses.adults*1.3,doses.adults*1.6)
dose.tox.weak.weaker <- function(doses){dose.tox.weak(doses)^1.3}; dose.tox.weak.stronger <- function(doses){dose.tox.weak(doses)^0.7}
df <- data.frame(doses.adults=rep(c(doses.adults),3),
                 doses.ped=rep(c(doses.ped),3),
                 y=c(dose.tox.weak(doses.ped),dose.tox.weak.weaker(doses.ped),dose.tox.weak.stronger(doses.ped)))
p1 <- ggplot(data.frame(x=rep(c(doses.adults),3),y=c(dose.tox.weak(doses.adults),dose.tox.weak.weaker(doses.adults),dose.tox.weak.stronger(doses.adults))),aes(x,y)) + 
  stat_function(fun = dose.tox.weak, aes(linetype = "Same as adults")) + 
  stat_function(fun = dose.tox.weak.weaker, aes(linetype = "Stronger")) + 
  stat_function(fun = dose.tox.weak.stronger, aes(linetype = "Weaker")) + 
  scale_linetype_manual(name="Pediatric toxicity scenario", values=c("solid", "dotted", "dashed")) +
  geom_point() + geom_point(data=df[!df$doses.ped==doses.adults,], aes(x=doses.ped, y=y), shape=1, size=2) + 
  scale_x_continuous(expand = c(0, 0.1), limits=c(0,80)) + scale_y_continuous(expand = c(0, 0.02), limits=c(0,1)) + 
  theme_bw() + ggtitle("2P logistic function, weak toxicity") +
  xlab("Dose") + ylab("Toxicity probability") + annotate("text", x=30, y=0.8, label="b = 1.5, e = 70, a = {0.7, 1, 1.3}", size=3)


dose.tox.moderate.weaker <- function(doses){dose.tox.moderate(doses)^1.3}; dose.tox.moderate.stronger <- function(doses){dose.tox.moderate(doses)^0.7}
df <- data.frame(doses.adults=rep(c(doses.adults),3),
                 doses.ped=rep(c(doses.ped),3),
                 y=c(dose.tox.moderate(doses.ped),dose.tox.moderate.weaker(doses.ped),dose.tox.moderate.stronger(doses.ped)))
p2 <- ggplot(data.frame(x=rep(c(doses.adults),3),y=c(dose.tox.moderate(doses.adults),dose.tox.moderate.weaker(doses.adults),dose.tox.moderate.stronger(doses.adults))),aes(x,y)) + 
  stat_function(fun = dose.tox.moderate, aes(linetype = "Same as adults")) + 
  stat_function(fun = dose.tox.moderate.weaker, aes(linetype = "Stronger")) + 
  stat_function(fun = dose.tox.moderate.stronger, aes(linetype = "Weaker")) + 
  scale_linetype_manual(name="Pediatric toxicity scenario", values=c("solid", "dotted", "dashed")) +
  geom_point() + geom_point(data=df[!df$doses.ped==doses.adults,], aes(x=doses.ped, y=y), shape=1, size=2) + scale_x_continuous(expand = c(0, 0.1), limits=c(0,80)) + scale_y_continuous(expand = c(0, 0.02), limits=c(0,1)) + 
  theme_bw() + ggtitle("2P logistic function, moderate toxicity") +
  xlab("Dose") + ylab("Toxicity probability") + annotate("text", x=30, y=0.8, label="b = 1.3, e = 50, a = {0.7, 1, 1.3}", size=3)

dose.tox.strong.weaker <- function(doses){dose.tox.strong(doses)^1.3}; dose.tox.strong.stronger <- function(doses){dose.tox.strong(doses)^0.4}
df <- data.frame(doses.adults=rep(c(doses.adults),3),
                 doses.ped=rep(c(doses.ped),3),
                 y=c(dose.tox.strong(doses.ped),dose.tox.strong.weaker(doses.ped),dose.tox.strong.stronger(doses.ped)))
p3 <- ggplot(data.frame(x=rep(c(doses.adults),3),y=c(dose.tox.strong(doses.adults),dose.tox.strong.weaker(doses.adults),dose.tox.strong.stronger(doses.adults))),aes(x,y)) + 
  stat_function(fun = dose.tox.strong, aes(linetype = "Same as adults")) + stat_function(fun = dose.tox.strong.weaker, aes(linetype = "Stronger")) + stat_function(fun = dose.tox.strong.stronger, aes(linetype = "Weaker")) + 
  scale_linetype_manual(name="Pediatric toxicity scenario", values=c("solid", "dotted", "dashed")) +
  geom_point() + geom_point(data=df[!df$doses.ped==doses.adults,], aes(x=doses.ped, y=y), shape=1, size=2) + scale_x_continuous(expand = c(0, 0.1), limits=c(0,80)) + scale_y_continuous(expand = c(0, 0.02), limits=c(0,1)) + 
  theme_bw() + ggtitle("2P logistic function, strong toxicity") +
  xlab("Dose") + ylab("Toxicity probability")  + annotate("text", x=45, y=0.1, label="b = 1.8, e = 27, a = {0.4, 1, 1.3}", size=3)



probs1 <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)
spline.plateau <- splinefun(c(3.5, doses.adults,188,40), probs1, "monoH.FC")
probs2 <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)^0.7
spline.plateau.stronger <- splinefun(c(3.5, doses.adults,188,40), probs2, "monoH.FC")
probs3 <- c(0.0001, 0.01, 0.03, 0.07, 0.10, 0.12, 0.15, 0.3, 0.95, 0.28)^1.3
spline.plateau.weaker <- splinefun(c(3.5, doses.adults,188,40), probs3, "monoH.FC")
df <- data.frame(doses.adults=rep(c(doses.adults),3),
                 doses.ped=rep(c(doses.ped),3),
                 y=c(spline.plateau(doses.ped),spline.plateau.weaker(doses.ped),spline.plateau.stronger(doses.ped)))
p4 <- ggplot(data.frame(x=rep(c(3.5, doses.adults,188),3),y=c(probs1[-9],probs2[-9],probs3[-9])),aes(x,y)) + 
  stat_function(fun = spline.plateau, aes(linetype = "Same as adults")) + stat_function(fun = spline.plateau.stronger, aes(linetype = "Stronger")) + stat_function(fun = spline.plateau.weaker, aes(linetype = "Weaker")) + 
  scale_linetype_manual(name="Pediatric toxicity scenario", values=c("solid", "dotted", "dashed")) +
  geom_point() + geom_point(data=df[!df$doses.ped==doses.adults,], aes(x=doses.ped, y=y), shape=1, size=2) + scale_x_continuous(expand = c(0, 0.1), limits=c(0,80)) + scale_y_continuous(expand = c(0, 0.02), limits=c(0,1)) +
  theme_bw() + ggtitle("Spline function, plateau form") +
  xlab("Dose") + ylab("Toxicity probability") + annotate("text", x=33, y=0.8, label="monotone Hermite spline, a = {0.7, 1, 1.3}", size=3)



probs1 <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)
spline.waves <- splinefun(c(doses.adults,188), probs1, "monoH.FC")
probs2 <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)^0.7
spline.waves.stronger <- splinefun(c(doses.adults,188), probs2, "monoH.FC")
probs3 <- c(0.01, 0.02, 0.13, 0.15, 0.27, 0.3, 0.5, 0.9)^1.3
spline.waves.weaker <- splinefun(c(doses.adults,188), probs3, "monoH.FC")
df <- data.frame(doses.adults=rep(c(doses.adults),3),
                 doses.ped=rep(c(doses.ped),3),
                 y=c(spline.waves(doses.ped),spline.waves.weaker(doses.ped),spline.waves.stronger(doses.ped)))

p5 <- ggplot(data.frame(x=rep(c(doses.adults,188),3),y=c(probs1,probs2,probs3)),aes(x,y)) + 
  stat_function(fun = spline.waves, aes(linetype = "Same as adults")) + stat_function(fun = spline.waves.stronger, aes(linetype = "Stronger")) + stat_function(fun = spline.waves.weaker, aes(linetype = "Weaker")) + 
  scale_linetype_manual(name="Pediatric toxicity scenario", values=c("solid", "dotted", "dashed")) +
  geom_point() + geom_point(data=df[!df$doses.ped==doses.adults,], aes(x=doses.ped, y=y), shape=1, size=2) + scale_x_continuous(expand = c(0, 0.1), limits=c(0,80)) + scale_y_continuous(expand = c(0, 0.02), limits=c(0,1)) +
  theme_bw() + ggtitle("Spline function, wave form") +
  xlab("Dose") + ylab("Toxicity probability") + annotate("text", x=33, y=0.8, label="monotone Hermite spline, a = {0.7, 1, 1.3}", size=3)



ggarrange(p1,p2,p3,p5,p4, common.legend=TRUE, legend = "bottom", labels = letters)

ggsave("scenarios.eps", device="eps", width=12, height=6, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")












################################################################################################################################################################################################################################################################################################################################################
# how many acceptable doses in the trials? 

nr.acceptable.trials <- function(tox.adults){
  newlabels <- c("same"="Same \nas adults",
                 "stronger"="Stronger \nthan adults",
                 "weaker"="Weaker \nthan adults")
  positions <- c("weaker", "same", "stronger")
  
  title <- if(tox.adults=="weak"){"\nAdult scenario:\nweakly toxic."} else if(tox.adults=="moderate"){"\nAdult scenario:\nmoderately toxic."} else if(tox.adults=="strong"){"\nAdult scenario:\nstrongly toxic."} else if(tox.adults=="spline.plateau"){"Adult scenario: \nspline function with\nlong plateau."} else if(tox.adults=="spline.waves"){"Adult scenario: \nspline function\nwith waves."}
  
  plotdata <- data %>% 
    filter(adult.scenario==tox.adults & !(specification %in% c("2P.weighted.25", "2P.weighted.50", "2P.weighted.75"))) %>%
    group_by(ped.scenario) %>% summarize(acceptable = sum(true.tox1<(1/3) & true.tox4>(1/6)), nonacceptable = sum((true.tox1<(1/3) & true.tox4>(1/6))==FALSE)) %>% 
    pivot_longer(., cols=c(acceptable, nonacceptable), names_to="acceptable", values_to="n") %>% group_by(ped.scenario) %>% mutate(percentage = n / sum(n)) %>% add_column("range"="all.doses")
  
  plotdata2 <- data %>% 
    filter(adult.scenario==tox.adults & !(specification %in% c("2P.weighted.25", "2P.weighted.50", "2P.weighted.75"))) %>%
    group_by(ped.scenario) %>% summarize(acceptable = sum(true.tox2<(1/3) & true.tox2>(1/6)), nonacceptable = sum((true.tox2<(1/3) & true.tox2>(1/6))==FALSE)) %>% 
    pivot_longer(., cols=c(acceptable, nonacceptable), names_to="acceptable", values_to="n") %>% group_by(ped.scenario) %>% mutate(percentage = n / sum(n)) %>% add_column("range"="dose2")
  
  p <- rbind(plotdata, plotdata2) %>% filter(acceptable=="acceptable") %>%
    ggplot(.) + geom_col(aes(x=ped.scenario, y=percentage, fill=range), position="dodge") + 
    scale_x_discrete(labels=newlabels, limits=positions) + 
    scale_fill_manual(name = "", labels = c("\nTrials with\nany acceptable dose\n", "\nTrials with\nacceptable dose 2\n"), values=c("grey70","grey30")) +  
    xlab("Pediatric toxicity scenario") + ylab("Proportion") + scale_y_continuous(limits=c(0,1)) + 
    theme_bw() + ggtitle(title)
  
  return(p)
}

ggarrange(
  nr.acceptable.trials(tox.adults = "weak"),
  nr.acceptable.trials(tox.adults = "moderate"),
  nr.acceptable.trials(tox.adults = "strong"),
  nr.acceptable.trials(tox.adults = "spline.waves"),
  nr.acceptable.trials(tox.adults = "spline.plateau"),
  ncol=5,nrow=1, common.legend=TRUE,legend="bottom")
ggsave("acceptable_trials.eps", device="eps", width=15, height=5, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")


data1 <- data %>% 
  group_by(adult.scenario, ped.scenario) %>% 
  summarize(acceptable = sum(true.tox1<(1/3) & true.tox4>(1/6)), nonacceptable = sum((true.tox1<(1/3) & true.tox4>(1/6))==FALSE)) %>% 
  pivot_longer(., cols=c(acceptable, nonacceptable), names_to="toxicity", values_to="n") %>% 
  group_by(adult.scenario, ped.scenario) %>% mutate(percentage = n / sum(n)) %>% add_column("column"="acceptable.all.doses") %>% filter(toxicity=="acceptable")

data2 <- data %>% 
  group_by(adult.scenario, ped.scenario) %>% 
  summarize(acceptable = sum(true.tox2<(1/3) & true.tox2>(1/6)), nonacceptable = sum((true.tox2<(1/3) & true.tox2>(1/6))==FALSE)) %>% 
  pivot_longer(., cols=c(acceptable, nonacceptable), names_to="toxicity", values_to="n") %>% 
  group_by(adult.scenario, ped.scenario) %>% mutate(percentage = n / sum(n)) %>% add_column("column"="acceptable.dose2") %>% filter(toxicity=="acceptable") 

data3 <- data %>% 
  group_by(adult.scenario, ped.scenario) %>% 
  summarize(not.too.toxic = sum(true.tox1<(1/3)), too.toxic = sum(true.tox1>(1/3))) %>% 
  pivot_longer(., cols=c(not.too.toxic, too.toxic), names_to="toxicity", values_to="n") %>% 
  group_by(adult.scenario, ped.scenario) %>% mutate(percentage = n / sum(n)) %>% add_column("column"="all.too.toxic") %>% filter(toxicity=="too.toxic") 

data4 <- data %>% 
  group_by(adult.scenario, ped.scenario) %>% 
  summarize(not.all.underdosing = sum(true.tox4>(1/6)), all.underdosing = sum(true.tox4<(1/6))) %>% 
  pivot_longer(., cols=c(all.underdosing, not.all.underdosing), names_to="toxicity", values_to="n") %>% 
  group_by(adult.scenario, ped.scenario) %>% mutate(percentage = n / sum(n)) %>% add_column("column"="all.underdosing") %>% filter(toxicity=="all.underdosing") 

rbind(data1, data2, data3, data4) %>% mutate(percentage = percentage*100) %>% pivot_wider(., id_cols=c(adult.scenario, ped.scenario), names_from=column, values_from=percentage) %>%
  xtable() %>% print(booktabs=TRUE)





################################################################################################################################################################################################################################################################################################################################################
# accuracy

newlabels <- c("1P.weak"="1P\nweak","1P.calibrated"="1P\ncalib","1P.full"="1P\nfull","1P.mixture"="1P\nmixt","1P.oracle"="1P\noracle",
               "2P.calibrated"="2P\ncalib","2P.full"="2P\nfull","2P.weak"="2P\nweak","2P.mixture"="2P\nmixt","2P.oracle"="2P\noracle")

positions.adult.scenario <- c("weak", "moderate", "strong", "spline.plateau", "spline.waves")
positions.ped.scenario <- c("weaker", "same", "stronger")

# see ?labeller
scenario_labeller <- labeller(
  `Pediatric toxicity` = c(`weaker` = "Pediatric toxicity:\nweaker", `same` = "Pediatric toxicity:\nsame", `stronger` = "Pediatric toxicity:\nstronger"),
  `Adult toxicity` = c(`weak` = "Adult toxicity:\nweak", `moderate` = "Adult toxicity:\nmoderate", `strong` = "Adult toxicity:\nstrong",
                       `spline.plateau` = "Adult toxicity:\nspline function with plateau form", `spline.waves` = "Adult toxicity:\nspline function with wave form"),
  .default = label_both
) 


plot1 <- data %>% 
  group_by("Prior" = unlist(specification), adult.scenario, ped.scenario) %>% 
  summarize("Acceptable" = sum(toxicity.MTD=="acceptable"),
            "Correct ET" = sum(toxicity.MTD=="correct ET"),
            "Sum of \ncorrect decisions" = sum(toxicity.MTD=="acceptable" | toxicity.MTD=="correct ET"),
            "N"=n()) %>% 
  pivot_longer(., cols=c(Acceptable, `Correct ET`, `Sum of \ncorrect decisions`), names_to="Decision", values_to="Number") %>%
  rowwise() %>%
  mutate("Proportion" = 100*Number/N, "lower" = prop.test(Number,N)$conf.int[1]*100, "upper" = prop.test(Number,N)$conf.int[2]*100,
         "Adult toxicity" = factor(adult.scenario, levels=positions.adult.scenario),
         "Pediatric toxicity" = factor(ped.scenario, levels=positions.ped.scenario)) %>% 
  ggplot(., aes(Prior, Proportion, color=Decision)) + geom_point(size=0.5, alpha = 1, position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0,100), breaks = seq(0,100,by=20), minor_breaks = seq(0,100,by=10)) + 
  xlab("Algorithm and prior configuration") + 
  scale_x_discrete(labels=newlabels, limits=positions, position=c("bottom")) + 
  scale_color_viridis_d(begin = 0.1, end = 0.9) +
  ylab("Percentage") + #ggtitle(title) +
  geom_errorbar(aes(Prior, y=Proportion, ymin=lower, ymax=upper), size=0.5, position=position_dodge(width=0.5)) + 
  geom_segment(aes(x = c(rep(0.5,15), rep(5.5,15)), y = Proportion, xend = c(rep(5.5,15), rep(10.5, 15)), yend = Proportion),
               data = . %>% filter(Decision == 'Sum of \ncorrect decisions', Prior %in% c('1P.oracle', '2P.oracle')),
               color="black", linetype="dashed", size=0.5) +
  facet_grid(`Pediatric toxicity` ~ `Adult toxicity`, labeller = scenario_labeller) +
  theme_bw() + theme(axis.text.x = element_text(size=7), legend.position = "bottom", panel.grid.major.x = element_blank()) +
  geom_vline(xintercept=seq(1.5,9.5), colour="grey90") # remove all grid lines on x axis and add new ones

library(gridExtra)
table <- t(c("'1P': 1P-CRM  ", "'2P': 2P-CRM  ", "'weak': weakly informative prior  ", "'mixt': mixture prior  ", "'calib': calibrated prior  ", "'full': full borrowing  "))
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 8,padding = unit(c(2, 4), "mm"))
tbl <- tableGrob(table, rows=NULL, theme=tt)

grid.arrange(plot1, tbl, nrow = 2, heights = c(4,0.2))
plot <- arrangeGrob(plot1, tbl, nrow = 2, heights = c(4,0.2))

ggsave("accuracy_all.eps", plot, device="eps", width=17, height=8, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")



################################################################################################################################################################################################################
# summary plots: if second dose was suitable

summary.plot <- function(tox.adults){
  
  title <- if(tox.adults=="weak"){"Adult scenario: weak toxicity."} else if(tox.adults=="moderate"){"Adult scenario: moderate toxicity."} else if(tox.adults=="strong"){"Adult scenario: strong toxicity."} else if(tox.adults=="spline.waves"){"Adult scenario: spline function with wave form."} else if(tox.adults=="spline.plateau"){"Adult scenario: spline function with long plateau."}
  p1 <- data %>% 
    filter(adult.scenario==tox.adults & ped.scenario=="same" & (skeleton2>0.2 & skeleton2<=0.3) & (tox.dose2>0.2 & tox.dose2<=0.3) & !(specification %in% c("2P.weighted.25", "2P.weighted.50", "2P.weighted.75"))) %>%
    group_by("Prior" = unlist(specification)) %>% 
    summarize("Acceptable" = sum(toxicity.MTD=="acceptable"),
              "Correct ET" = sum(toxicity.MTD=="correct ET"),
              "Sum of \ncorrect decisions" = sum(toxicity.MTD=="acceptable" | toxicity.MTD=="correct ET"),
              "N"=n()) %>% 
    pivot_longer(., cols=c(Acceptable, `Correct ET`, `Sum of \ncorrect decisions`), names_to="Decision", values_to="Number") %>%
    rowwise() %>%
    mutate("Proportion" = 100*Number/N, "lower" = prop.test(Number,N)$conf.int[1]*100, "upper" = prop.test(Number,N)$conf.int[2]*100) %>%
    ggplot(., aes(Prior, Proportion, color=Decision)) + geom_point(size=0.5, alpha = 1, position=position_dodge(width=0.5)) +
    scale_y_continuous(limits=c(0,100)) + xlab("Algorithm and prior configuration") + scale_x_discrete(labels=newlabels, limits=positions) + 
    scale_color_viridis_d() +
    ylab("Percentage of\nacceptable MTDs") + ggtitle(title) +
    geom_errorbar(aes(Prior, y=Proportion, ymin=lower, ymax=upper), size=0.5, position=position_dodge(width=0.5)) +
    geom_segment(aes(x = c(0.5, 5.5), y = Proportion, xend = c(5.5, 10.5), yend = Proportion), 
                 data = . %>% filter(Decision == 'Sum of \ncorrect decisions', Prior %in% c('1P.oracle', '2P.oracle')),
                 color="black", linetype="dashed", size=0.5) + 
    theme_bw()
  return(p1)
}
summary.plot(tox.adults="weak")
# The second pediatric dose is an acceptable dose, i.e. the borrowed information on the MTD was suitable.
p1 <- summary.plot(tox.adults="weak")
p2 <- summary.plot(tox.adults="moderate")
p3 <- summary.plot(tox.adults="strong")
p4 <- summary.plot(tox.adults="spline.waves")
p5 <- summary.plot(tox.adults="spline.plateau")
ggarrange(p1, p2, p3, p4, p5, ncol=1, nrow=5)
ggsave("accuracy_2nd_dose.eps", device="eps", width=8, height=10, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")


########################################################################################################################################################################################################
# Number of DLTs
positions <- c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full")

data$DLTs <- data$toxicities.per.dose1 + data$toxicities.per.dose2 + data$toxicities.per.dose3 + data$toxicities.per.dose4

data$specification <- unlist(data$specification)

prior_names <- c(
  `1P.weak` = "1P-CRM with \nweak prior", `1P.mixture` = "1P-CRM with \nmixture prior", `1P.calibrated` = "1P-CRM with \ncalibrated prior",
  `1P.full` = "1P-CRM with \nfull borrowing", `1P.oracle` = "1P-CRM \noracle", `2P.weak` = "2P-CRM with\nweak prior", `2P.mixture` = "2P-CRM with \nmixture prior",
  `2P.calibrated` = "2P-CRM with \ncalibrated prior", `2P.full` = "2P-CRM with \nfull borrowing", `2P.oracle` = "2P-CRM \noracle")


scenario_labeller <- labeller(
  ped.scenario = c(`weaker` = "Pediatric toxicity:\nweaker", `same` = "Pediatric toxicity:\nsame", `stronger` = "Pediatric toxicity:\nstronger"),
  specification = c(
    `1P.weak` = "1P-CRM with \nweak prior", `1P.mixture` = "1P-CRM with \nmixture prior", `1P.calibrated` = "1P-CRM with \ncalibrated prior",
    `1P.full` = "1P-CRM with \nfull borrowing", `1P.oracle` = "1P-CRM \noracle", `2P.weak` = "2P-CRM with\nweak prior", `2P.mixture` = "2P-CRM with \nmixture prior",
    `2P.calibrated` = "2P-CRM with \ncalibrated prior", `2P.full` = "2P-CRM with \nfull borrowing", `2P.oracle` = "2P-CRM \noracle"
  ),
  .default = label_both
) 

positions.ped.scenario <- c("weaker", "same", "stronger")

data %>%
  filter(adult.scenario=="moderate" & sample.size==10 &
           (specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full"))) %>%
  group_by(ped.scenario, specification) %>% 
  summarise("0" = sum(DLTs == 0),"1" = sum(DLTs == 1),"2" = sum(DLTs == 2),"3" = sum(DLTs == 3),
            "4" = sum(DLTs == 4),"5" = sum(DLTs == 5),"6" = sum(DLTs == 6),"7" = sum(DLTs == 7),
            "8" = sum(DLTs == 8)) %>% 
  gather(DLTs, n, -c(specification, ped.scenario)) %>% 
  mutate(specification = factor(specification, levels=positions),
         ped.scenario = factor(ped.scenario, levels=positions.ped.scenario)) %>% 
  group_by(ped.scenario, specification) %>% 
  mutate(groupmean = sum(as.numeric(DLTs)*n)/sum(n), sum_n=sum(n), 
         prop = 100 * n / sum(n), 
         DLTs = factor(DLTs, levels=as.character(0:10))) %>% #print(n=nrow(.))
  ggplot(., aes(x = DLTs, y = prop, fill = DLTs)) +
  geom_bar(stat="identity", position = "dodge", width=0.7) + 
  scale_y_continuous(limits=c(0,50), breaks=seq(0,50,by=20), minor_breaks = seq(0,50,by=10)) + 
  xlab("Number of DLTs") + ylab("Percentage of simulated pediatric trials") +
  coord_flip() + facet_grid(ped.scenario ~ specification, labeller = scenario_labeller) +
  theme_bw() + geom_vline(aes(xintercept=groupmean+1), linetype="dashed", size=1) +
  theme(legend.position = "none",strip.text.x = element_text(size = 12),axis.title = element_text(size=12)) + 
  guides(colour=FALSE) + scale_fill_viridis_d(end=0.9)
ggsave("DLTs.eps", device="eps", width=14, height=8,dpi=800, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")




########################################################################################################################################################################################################



# patient allocation I: patients at acceptably toxic doses
patients.acceptable1 <- ifelse(data$true.tox1<(1/3) & data$true.tox1>(1/6), data$patients.per.dose1, 0)
patients.acceptable2 <- ifelse(data$true.tox2<(1/3) & data$true.tox2>(1/6), data$patients.per.dose2, 0)
patients.acceptable3 <- ifelse(data$true.tox3<(1/3) & data$true.tox3>(1/6), data$patients.per.dose3, 0)
patients.acceptable4 <- ifelse(data$true.tox4<(1/3) & data$true.tox4>(1/6), data$patients.per.dose4, 0)
data$patients.acceptable <- patients.acceptable1+patients.acceptable2+patients.acceptable3+patients.acceptable4

# was there any acceptable dose? otherwise we cannot allocate patients to acceptable doses anyway
data$any.acceptable.dose <- ifelse((data$true.tox1<(1/3) & data$true.tox1>(1/6)) | data$true.tox2<(1/3) & data$true.tox2>(1/6) | data$true.tox3<(1/3) & data$true.tox3>(1/6) | data$true.tox4<(1/3) & data$true.tox4>(1/6), 
                                   1, 0)

scenario_labeller <- labeller(
  ped.scenario = c(`weaker` = "Pediatric toxicity:\nweaker", `same` = "Pediatric toxicity:\nsame", `stronger` = "Pediatric toxicity:\nstronger"),
  specification = c(
    `1P.weak` = "1P-CRM with \nweak prior", `1P.mixture` = "1P-CRM with \nmixture prior", `1P.calibrated` = "1P-CRM with \ncalibrated prior",
    `1P.full` = "1P-CRM with \nfull borrowing", `1P.oracle` = "1P-CRM \noracle", `2P.weak` = "2P-CRM with\nweak prior", `2P.mixture` = "2P-CRM with \nmixture prior",
    `2P.calibrated` = "2P-CRM with \ncalibrated prior", `2P.full` = "2P-CRM with \nfull borrowing", `2P.oracle` = "2P-CRM \noracle"
  ),
  .default = label_both
) 

positions.ped.scenario <- c("weaker", "same", "stronger")

data %>%
  filter(adult.scenario=="moderate" & sample.size==10 & any.acceptable.dose==1 &
           (specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full"))) %>%
  group_by(ped.scenario, specification) %>% 
  summarise("0" = sum(patients.acceptable == 0),
            "2" = sum(patients.acceptable == 2),
            "4" = sum(patients.acceptable == 4),
            "6" = sum(patients.acceptable == 6),
            "8" = sum(patients.acceptable == 8),
            "10" = sum(patients.acceptable == 10)) %>% 
  gather(patients.acceptable, n, -c(specification, ped.scenario)) %>% 
  mutate(specification = factor(specification, levels=positions),
         ped.scenario = factor(ped.scenario, levels=positions.ped.scenario)) %>% 
  group_by(ped.scenario, specification) %>% 
  mutate(groupmean = 0.5 * sum(as.numeric(patients.acceptable)*n)/sum(n), sum_n=sum(n), # groupmean * 0.5, because the discrete scale is in stepts of 2
         prop = 100 * n / sum(n), 
         patients.acceptable = factor(patients.acceptable, levels=c('0', '2', '4', '6', '8', '10'))) %>% #print(n=nrow(.))
  ggplot(., aes(x = patients.acceptable, y = prop, fill = patients.acceptable)) +
  geom_bar(stat="identity", position = "dodge", width=0.7) + 
  scale_y_continuous(limits=c(0,65), breaks=seq(0,65,by=20), minor_breaks = seq(0,65,by=10)) + 
  xlab("Number of Patients\nat acceptable doses") + ylab("Percentage of simulated pediatric trials with at least one acceptable dose") +
  coord_flip() + facet_grid(ped.scenario ~ specification, labeller = scenario_labeller) +
  theme_bw() + geom_vline(aes(xintercept=groupmean+1), linetype="dashed", size=1) +
  theme(legend.position = "none",strip.text.x = element_text(size = 12),axis.title = element_text(size=12)) + 
  guides(colour=FALSE) + scale_fill_viridis_d(end=0.9)
ggsave("allocation.eps", device="eps", width=14, height=8, dpi=800, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")


# get raw numbers for text
tox.adults="moderate"; tox.ped="weaker"
data %>%
  filter(adult.scenario==tox.adults & ped.scenario==tox.ped & sample.size==10 & any.acceptable.dose==1 &
           (specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full"))) %>%
  group_by(specification) %>% 
  summarise("0" = sum(patients.acceptable == 0),
            "2" = sum(patients.acceptable == 2),
            "4" = sum(patients.acceptable == 4),
            "6" = sum(patients.acceptable == 6),
            "8" = sum(patients.acceptable == 8),
            "10" = sum(patients.acceptable == 10),
            n_total=n()) %>% 
  pivot_longer(data=., cols=c('0', '2', '4', '6', '8', '10'), names_to="Patients") %>% mutate(prop=value/n_total*100) %>% 
  filter(Patients==0) %>% print(n=nrow(.))
  


###################################################################################################################
# Oscillation: are there trials which recommend a dose on which no patient was treated?
data$oscillation <- ifelse(data$MTD.dose.level==0, FALSE,
                           ifelse(data$MTD.dose.level==1, data$patients.per.dose1==0,
                                  ifelse(data$MTD.dose.level==2, data$patients.per.dose2==0,
                                         ifelse(data$MTD.dose.level==3, data$patients.per.dose3==0,
                                                ifelse(data$MTD.dose.level==4, data$patients.per.dose4==0,NA)))))

positions <- c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full")



tox.adults <- "moderate"
tox.ped <- "same"

p1 <- data %>% 
  filter(adult.scenario==tox.adults & ped.scenario==tox.ped & sample.size==10) %>% 
  filter(specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full")) %>% 
  group_by("Prior" = specification) %>% 
  summarize("Oscillation" = sum(oscillation)/n()*100) %>% 
  ggplot(., aes(Prior, Oscillation, color=Prior)) + geom_point(size=4, stroke=2) +
  scale_shape_manual(values=seq(0,9)) +
  scale_y_continuous(limits=c(0,15)) + ylab("Percentage of trials with oscillation") +
  scale_x_discrete(labels=newlabels, limits=positions) + xlab("Algorithm and prior configuration") + 
  ggtitle("Adult toxicity moderate. \nPediatric toxicity the same.") +
  theme_bw() + theme(legend.position = "none",plot.title = element_text(size=10),axis.title = element_text(size=10)) + scale_color_viridis_d()


tox.adults <- "moderate"
tox.ped <- "weaker"

p2 <- data %>% 
  filter(adult.scenario==tox.adults & ped.scenario==tox.ped & sample.size==10) %>% 
  filter(specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full")) %>% 
  group_by("Prior" = specification) %>% 
  summarize("Oscillation" = sum(oscillation)/n()*100) %>% 
  ggplot(., aes(Prior, Oscillation, color=Prior)) + geom_point(size=4, stroke=2)  +
  scale_shape_manual(values=seq(0,9)) +
  scale_y_continuous(limits=c(0,15)) + ylab("Percentage of trials with oscillation") +
  scale_x_discrete(labels=newlabels, limits=positions) + xlab("Algorithm and prior configuration") + 
  ggtitle("Adult toxicity moderate. \nPediatric toxicity weaker.") +
  theme_bw() + theme(legend.position = "none",plot.title = element_text(size=10),axis.title.y = element_blank()) + scale_color_viridis_d()



tox.adults <- "strong"
tox.ped <- "same"

p3 <- data %>% 
  filter(adult.scenario==tox.adults & ped.scenario==tox.ped & sample.size==10) %>% 
  filter(specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full")) %>% 
  group_by("Prior" = specification) %>% 
  summarize("Oscillation" = sum(oscillation)/n()*100) %>% 
  ggplot(., aes(Prior, Oscillation, color=Prior)) + geom_point(size=4, stroke=2) +
  scale_shape_manual(values=seq(0,9)) +
  scale_y_continuous(limits=c(0,15)) + ylab("Percentage of trials with oscillation") +
  scale_x_discrete(labels=newlabels, limits=positions) + xlab("Algorithm and prior configuration") + 
  ggtitle("Adult toxicity strong \nPediatric toxicity the same.") +
  theme_bw() + theme(legend.position = "none",plot.title = element_text(size=10),axis.title = element_text(size=10)) + scale_color_viridis_d()


tox.adults <- "strong"
tox.ped <- "weaker"

p4 <- data %>% 
  filter(adult.scenario==tox.adults & ped.scenario==tox.ped & sample.size==10) %>% 
  filter(specification %in% c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full")) %>% 
  group_by("Prior" = specification) %>% 
  summarize("Oscillation" = sum(oscillation)/n()*100) %>% 
  ggplot(., aes(Prior, Oscillation, color=Prior)) + geom_point(size=4, stroke=2)  +
  scale_shape_manual(values=seq(0,9)) +
  scale_y_continuous(limits=c(0,15)) + ylab("Percentage of trials with oscillation") +
  scale_x_discrete(labels=newlabels, limits=positions) + xlab("Algorithm and prior configuration") + 
  ggtitle("Adult toxicity strong. \nPediatric toxicity weaker.") +
  theme_bw() + theme(legend.position = "none",plot.title = element_text(size=10),axis.title.y = element_blank()) + scale_color_viridis_d()


ggarrange(p1,p2,p3,p4,nrow=1)
ggsave("oscillation.eps", device="eps", width=16, height=6, path = "S:/C01/iBikE/Studien/Phase1Studien/3_Programme/R-Codes/Figures")




###############################################################################################################################################################################################################
# summary table



data %>% 
  filter(!(specification %in% c("2P.weighted.25", "2P.weighted.50", "2P.weighted.75"))) %>%
  group_by(adult.scenario, ped.scenario, "Prior" = specification) %>% 
  summarize("% with acceptable dose" = round(sum(any.acceptable.dose)/n()*100,1),
            "% correct decision" = round(sum(toxicity.MTD=="acceptable" | toxicity.MTD=="correct ET")/n()*100,1),
            "mean number of DLTs" = mean(DLTs), 
            "mean number of patients at acceptable doses" = mean(patients.acceptable),
            "% oscillation" = sum(oscillation)/n()*100) %>% 
  mutate(adult.scenario = factor(adult.scenario, levels=c("weak", "moderate", "strong", "spline.plateau", "spline.waves")),
         ped.scenario = factor(ped.scenario, levels=c("weaker", "same", "stronger")),
         Prior = factor(Prior, levels=c("1P.oracle","1P.weak","1P.mixture","1P.calibrated","1P.full","2P.oracle","2P.weak","2P.mixture","2P.calibrated","2P.full"))) %>% 
  arrange(adult.scenario, ped.scenario, Prior) %>% relocate(4, .after=2) %>%
  xtable(digits = c(0,0,0,1,1,1,2,2,1)) %>% print(booktabs=TRUE)




# The other exception was the spline function with wave form, when pediatric toxicity was weaker than in adults. In this case, all configurations achieved accuracy rates above 70 \%. This was accountable to an artefact of the simulation scenarios ... 
data %>% filter(adult.scenario=="spline.waves" & ped.scenario=="weaker" & specification=="1P.full")
data %>% filter(adult.scenario=="spline.waves" & ped.scenario=="weaker" & specification=="1P.full") %>% group_by(doses2) %>% count()
spline.waves.weaker(c(36,47))














