# Compare different model using ASPREE data
# effect of diabetes on the risk of physical disability
# Author: Thao Le

rm(list = ls())
set.seed(7)
library(survival)
library(SmoothHazard)
library(tidyverse)
library(icenReg)


# prepare data
load("aspree.Rdata")

# all
all.d <- aspree.d %>% 
  mutate(mp = (ttevent + lag.days)/2,
         obs.d = ifelse(!is.na(Withdrawn_DSR) & !is.na(Death_DSR), Death*as.numeric(Death_DSR<Withdrawn_DSR),Death),
         obs.dt = ifelse(obs.d==1, Death_DSR, LastPhysDisScreen_DSR),
         tt.cens.death = ifelse(event == 1, ttevent, obs.dt),
         tt.cens.lastvisit = ttevent,
         tt.cens.lastsecondvisit = ifelse(event == 1, ttevent, lag.days),
         tt.cens.lastnegativetest = ifelse(event == 1, ttevent,
                                           ifelse(max.trigger == 1, lag.days, ttevent)),
         tt.cens.death.mp = ifelse(event == 1, mp, obs.dt),
         tt.cens.lastvisit.mp = ifelse(event == 1, mp, ttevent),
         tt.cens.lastsecondvisit.mp = ifelse(event == 1, mp, lag.days),
         tt.cens.lastnegativetest.mp = ifelse(event == 1, mp,
                                              ifelse(max.trigger == 1, lag.days, ttevent)),
         
         # ic model
         ic.l.lastvisit = ifelse(event == 1, lag.days, ttevent),
         ic.l.last2visit = ifelse(event == 1, lag.days, lag.days),
         ic.l.lastnegvisit = ifelse(event == 1, lag.days,
                                    ifelse(max.trigger == 1, lag.days, ttevent)),
         ic.l.death = ifelse(event == 1, lag.days, obs.dt),
         
         ic.u = ifelse(event == 1, ttevent, Inf),
         
         # idm models
         
         L.LV = ifelse(event == 1, lag.days, ttevent),
         R.LV = ttevent,
         
         L.LSecV = lag.days,
         R.LSecV = ifelse(event == 1, ttevent, lag.days),
         
         L.LNegV = ifelse(event == 1, lag.days, 
                          ifelse(max.trigger == 1, lag.days, ttevent)),
         R.LNegV = ifelse(event == 1, ttevent,
                          ifelse(max.trigger == 1, lag.days, ttevent))
         
  )

# Fit model 0-->1
lastvisit <- coxph(Surv(tt.cens.lastvisit, event) ~ Diab_deriv, data = all.d)
lastsecondvisit <- coxph(Surv(tt.cens.lastsecondvisit, event) ~ Diab_deriv, data = all.d)
lasttest <- coxph(Surv(tt.cens.lastnegativetest, event) ~ Diab_deriv, data = all.d)
censordead <- coxph(Surv(tt.cens.death, event) ~ Diab_deriv, data = all.d)

lastvisit.mp <- coxph(Surv(tt.cens.lastvisit.mp, event) ~ Diab_deriv, data = all.d)
lastsecondvisit.mp <- coxph(Surv(tt.cens.lastsecondvisit.mp, event) ~ Diab_deriv, data = all.d)
lasttest.mp <- coxph(Surv(tt.cens.lastnegativetest.mp, event) ~ Diab_deriv, data = all.d)
censordead.mp <- coxph(Surv(tt.cens.death.mp, event) ~ Diab_deriv, data = all.d)

ic.lastvisit <- ic_sp(Surv(ic.l.lastvisit, ic.u, type = "interval2") ~ Diab_deriv, data = all.d,bs_samples = 500)
ic.last2visit <- ic_sp(Surv(ic.l.last2visit, ic.u, type = "interval2") ~ Diab_deriv, data = all.d, bs_samples = 500)
ic.lastnegvisit <- ic_sp(Surv(ic.l.lastnegvisit, ic.u, type = "interval2") ~ Diab_deriv, data = all.d, bs_samples = 500)
ic.death <- ic_sp(Surv(ic.l.death, ic.u, type = "interval2") ~ Diab_deriv, data = all.d, bs_samples = 500)


idmfit <- idm(formula01 = Hist(time = list(lag.days, ttevent), event = event) ~ Diab_deriv, 
              formula02 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
              formula12 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
              data = as.data.frame(all.d))

idmfit.LV <- idm(formula01 = Hist(time = list(L.LV, R.LV), event = event) ~ Diab_deriv, 
                 formula02 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
                 formula12 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
                 data = as.data.frame(all.d))
idmfit.LSecV <- idm(formula01 = Hist(time = list(L.LSecV, R.LSecV), event = event) ~ Diab_deriv, 
                    formula02 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
                    formula12 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
                    data = as.data.frame(all.d))

idmfit.LNegV <- idm(formula01 = Hist(time = list(L.LNegV, R.LNegV), event = event) ~ Diab_deriv, 
                    formula02 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
                    formula12 = Hist(time = obs.dt, event = obs.d) ~ Diab_deriv,
                    data = as.data.frame(all.d), maxiter = 500, CV = 1)

### transistion from 0 --> 2, competing risks
all.d <- all.d %>% 
  mutate(Death_DSR = ifelse(!is.na(Death_DSR), Death_DSR, LastDeathScreen_DSR),
         event.d = ifelse(event == 1, 0, Death),
         ttevent.mp = ifelse(event == 1, mp, Death_DSR),
         ttevent.fp = ifelse(event == 1, ttevent, Death_DSR))

mp.Cox <- coxph(Surv(ttevent.mp, event.d) ~ Diab_deriv, data = all.d)
fp.Cox <- coxph(Surv(ttevent.fp, event.d) ~ Diab_deriv, data = all.d)

save(lastvisit, lastsecondvisit, lasttest, censordead, lastvisit.mp, lastsecondvisit.mp, lasttest.mp, censordead.mp,
     idmfit, idmfit.LV, idmfit.LSecV, idmfit.LNegV, 
     ic.last2visit, ic.lastnegvisit, ic.lastvisit, ic.death,
     mp.Cox, fp.Cox,
     file = "./Interdata/2023-01-18_aspree-idm.Rdata")


