# Aim: Summary simulation results
# author: Thao Le

library(tidyverse)
library(gt)
library(gtsummary)
library(here)


####################
# length = 0.5
####################
load("./Interdata/2023-02-17_simresults-All-Length05.Rdata")

scenarios <- scenarios %>% mutate(scenario = row_number())
paraX05 <- paraX %>% 
  filter(!(Estimate == 0| abs(Estimate) >100 | is.na(Estimate) | is.na(SE))) %>% 
  mutate(scenario = as.numeric(scenario))
names(paraX05) <- c("scenario", "Model", "beta01_est", "se_est")

l05.X <- paraX05 %>% 
  left_join(scenarios, by = "scenario") %>% 
  mutate(bias = beta01_est - beta01,
         CIcontain = as.numeric(beta01_est - 1.96*se_est <= beta01 & beta01 <= beta01_est + 1.96*se_est),
         length = 0.5)


tab2 <- l05.X %>% 
  group_by(scenario, Model) %>% 
  summarise(
    n.sim = n(),
    mean.bias = mean(bias),
    mean.estimate = mean(beta01_est),
    MCse.bias = sqrt(sum((beta01_est - mean.estimate)^2)/(n.sim*(n.sim - 1))),
    sd.bias = sd(bias),
    mean.coverage = mean(CIcontain),
    empse = sqrt(sum((beta01_est - mean.estimate)^2)/(n.sim - 1)),
  ) %>% 
  left_join(scenarios, by = "scenario") %>% 
  mutate(upper.ci = mean.bias + 1.96*MCse.bias,
         lower.ci = mean.bias - 1.96*MCse.bias)

parW05 <- paraW %>% 
  mutate(scenario = as.numeric(scenario),
         length = 0.5) %>% 
  group_by(scenario, names) %>% 
  summarise(
    mean.lambda = mean(lambda),
    mean.gamma = mean(gamma)
  ) %>% 
  left_join(scenarios, by = "scenario")

## nevents
event05 <- events %>% 
  mutate(scenario = as.numeric(scenario)) %>% 
  left_join(scenarios, by = "scenario")

## n0
n005 <- n0 %>% 
  rename(Model = names) %>% 
  mutate(scenario = as.numeric(scenario)) %>% 
  left_join(scenarios, by = "scenario")

####################
# length = 1
####################
load("./Interdata/2023-02-17_simresults-All-Length10.Rdata")

scenarios <- scenarios %>% mutate(scenario = row_number())
paraX1 <- paraX %>%
  filter(!(Estimate == 0| abs(Estimate) >100 | is.na(Estimate) | is.na(SE))) %>%
  mutate(scenario = as.numeric(scenario))
names(paraX1) <- c("scenario", "Model", "beta01_est", "se_est")

l10.X <- paraX1 %>%
  left_join(scenarios %>% mutate(scenario = row_number()), by = "scenario") %>%
  mutate(bias = beta01_est - beta01,
         CIcontain = as.numeric(beta01_est - 1.96*se_est <= beta01 & beta01 <= beta01_est + 1.96*se_est),
         length = 1)
tab3 <- l10.X %>%
  group_by(scenario, Model) %>%
  summarise(
    n.sim = n(),
    mean.bias = mean(bias),
    mean.estimate = mean(beta01_est),
    MCse.bias = sqrt(sum((beta01_est - mean.estimate)^2)/(n.sim*(n.sim - 1))),
    sd.bias = sd(bias),
    mean.coverage = mean(CIcontain),
    empse = sqrt(sum((beta01_est - mean.estimate)^2)/(n.sim - 1)),
  ) %>% 
  left_join(scenarios, by = "scenario") %>% 
  mutate(upper.ci = mean.bias + 1.96*MCse.bias,
         lower.ci = mean.bias - 1.96*MCse.bias)
## nevents
event10 <- events %>%
  mutate(scenario = as.numeric(scenario)) %>%
  left_join(scenarios, by = "scenario")

## n0
n010 <- n0 %>%
  rename(Model = names) %>%
  mutate(scenario = as.numeric(scenario)) %>%
  left_join(scenarios, by = "scenario")

##################
# combine
##################
# save scenarios
save(scenarios, file = ("./Interdata/2023-02-20_Simulation-Senarios.Rdata"))

# combine all X
all.bias <- rbind(tab2, tab3)
save(all.bias, file = ("./Interdata/2023-02-20_Bias-Results.Rdata"))

# combine all events
all.events <- rbind(event05,
                    event10)
save(all.events, file = "./Interdata/2023-02-20_Nevents-Results.Rdata")




