# Run simulation (clock forward time)
# Scenarios: Visit length = 0.5
# author: Thao Le

# library
set.seed(05052022)

library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(SmoothHazard)
library(furrr)
library(mstate)
library(SurvRegCensCov)
library(icenReg)
source("01_Functions.R")
source("02_RunFuns.R")

# scenarios
# length = 0.5, 20 inspections

allscenarios <- crossing(Niter = 1000,
                         n = 1000,
                         beta01 = 0.5, 
                         beta02 = c(0, 0.5), 
                         beta12 = c(0, 0.5, -0.5),
                         shape01 = c(1, 1.5),
                         shape02 = c(1, 1.2, 1.5),
                         shape12 = c(1, 1.2, 1.5),
                         scale01 = 0.1,
                         scale02 = 0.1,
                         scale12 = 0.1,
                         cova.prob = 0.5, 
                         inspectLength = 0.5, 
                         inspections = 10,
                         false.positive.rate = c(0, 1),
                         min.cens = 1,
                         max.cens = 1000,
                         max.int = 500
)


scenarios <- allscenarios %>%
  filter(!(beta02 == 0.5 & beta12 == -0.5)) %>% 
  filter((shape01 == 1 & shape02 == 1 & shape12 == 1)|
           (shape01 == 1 & shape02 == 1.2 & shape12 == 1.5)|
        ( shape01 == 1.5 & shape02 == 1.2 & shape12 == 1.2)|
         (shape01 == 1.5 & shape02 == 1.5 & shape12 == 1.5))

# run simulation
plan(multisession)
system.time(tmp <- future_pmap(.l = scenarios, 
                               .f = simulation_nrun,
                               .options = furrr_options(seed = 1)))
# combine results
paraX <- bind_rows(lapply(tmp, function(x) x$paraX), .id = "scenario")
rownames(paraX) <- NULL
paraW <- bind_rows(lapply(tmp, function(x) x$paraW), .id = "scenario")
rownames(paraW) <- NULL
n0 <- bind_rows(lapply(tmp, function(x) x$n0), .id = "scenario")
rownames(n0) <- NULL
events <- bind_rows(lapply(tmp, function(x) x$events), .id = "scenario")

save(paraX, paraW, n0, events, scenarios, file = "2023-02-17_simresults-All-Length05.Rdata")

