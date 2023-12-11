# run functions 
# author: thao le 
# date: 1 Sep 2022
# Update: clock forward 


###################################### 
# Function to do 1 simulation run
# output: 
# - estimate for beta01
# - estimate of shape and scale for weibull distributions
# - number of events
###################################### 
simulation_1run <- function(n,
                            beta01, 
                            beta02, 
                            beta12,
                            shape01,
                            shape02,
                            shape12,
                            scale01,
                            scale02,
                            scale12,
                            cova.prob, 
                            inspectLength, 
                            inspections,
                            false.positive.rate,
                            min.cens = 0,
                            max.cens = 10,
                            max.int = 500
){
  d <- sim.illdeath.Interval.FP(n = n, 
                                shape01 = shape01, shape02 = shape02, shape12 = shape12,
                                scale01 = scale01, scale02 = scale02, scale12 = scale12, 
                                beta01 = beta01, beta02 = beta02, beta12 = beta12, 
                                cova.prob = cova.prob, inspections = inspections, 
                                inspectLength = inspectLength,
                                false.positive.rate =  false.positive.rate, 
                                min.cens = min.cens, max.cens = max.cens, 
                                max.int = max.int)
  # data analysis
  true.01 <- WeibullReg(Surv(illt, ills) ~ X,  data = d)
  
  # prepare event time for model fitting
  d <- d %>% 
    mutate(
      midpoint = (tstart + tstop)/2,
      
      # cox model
      cox.lastvisit = ifelse(persist.endpoint == 1, midpoint, tstop),
      cox.last2visit = ifelse(persist.endpoint == 1, midpoint, tstart),
      cox.lastnegvisit = ifelse(persist.endpoint == 1, midpoint,
                                ifelse(test.result == 1, tstart, tstop)),
      cox.deathtime = ifelse(persist.endpoint == 1, midpoint, 
                             ifelse(obs.ds == 1, obs.dt, tstop)),
      
      # Weibull ic
      ic.l.lastvisit = ifelse(persist.endpoint == 1, tstart, tstop),
      ic.l.last2visit = tstart,
      ic.l.lastnegvisit = ifelse(persist.endpoint == 1, tstart, 
                                 ifelse(test.result == 1, tstart, tstop)),
      ic.l.deathtime = ifelse(persist.endpoint == 1, tstart, 
                              ifelse(obs.ds == 1, obs.dt, tstop)),
      ic.u = ifelse(persist.endpoint == 1, tstop, Inf),
      
      # idm models: no event, l = r
      idm.l.lastvisit = ifelse(persist.endpoint == 1, tstart, tstop),
      idm.r.lastvisit = tstop,
      
      idm.l.last2visit = tstart,
      idm.r.last2visit = ifelse(persist.endpoint == 1, tstop, tstart),
      
      idm.l.lastnegvisit = ifelse(persist.endpoint == 1, tstart, 
                                  ifelse(test.result == 1, tstart, tstop)),
      idm.r.lastnegvisit = ifelse(persist.endpoint == 1, tstop, 
                                  ifelse(test.result == 1, tstart, tstop))
      
      
    )
  # survreg model
  survreg.fit <- function(d, time.var){
    newd <- d %>% filter(.data[[time.var]] >0)
    n0 = sum(d[[time.var]]==0)
    fit <- WeibullReg(Surv(newd[[time.var]], persist.endpoint) ~ X, 
                      data = newd)
    return(list(fit = fit, n0 = n0))
  }
  
  # model fit:
  # - Weibull model with mid point for event time
  W.lastvisit <- survreg.fit(d = d, time.var = "cox.lastvisit")
  W.last2visit <- survreg.fit(d = d, time.var = "cox.last2visit")
  W.lastnegvisit <- survreg.fit(d = d, time.var = "cox.lastnegvisit")
  W.deathtime <- survreg.fit(d = d, time.var = "cox.deathtime")
  
  # - IC Weibull model 
  ic.lastvisit <- ic_par(Surv(ic.l.lastvisit, ic.u, type = "interval2") ~ X, data = d)
  ic.last2visit <- ic_par(Surv(ic.l.last2visit, ic.u, type = "interval2") ~ X, data = d)
  ic.lastnegvisit <- ic_par(Surv(ic.l.lastnegvisit, ic.u, type = "interval2") ~ X, data = d)
  ic.deathtime <- ic_par(Surv(ic.l.deathtime, ic.u, type = "interval2") ~ X, data = d)
  
  # - IDM models 
  idm.lastvisit <- idm(formula01 = Hist(time = list(idm.l.lastvisit, idm.r.lastvisit), 
                                        event = persist.endpoint) ~ X, 
                       formula02 = Hist(time = obs.dt, event = obs.ds) ~ X,
                       formula12 = Hist(time = obs.dt, event = obs.ds) ~ X,
                       data = as.data.frame(d), maxiter = 500, CV = 1)
  idm.last2visit <- idm(formula01 = Hist(time = list(idm.l.last2visit, idm.r.last2visit), 
                                         event = persist.endpoint) ~ X, 
                        formula02 = Hist(time = obs.dt, event = obs.ds) ~ X,
                        formula12 = Hist(time = obs.dt, event = obs.ds) ~ X,
                        data = as.data.frame(d), maxiter = 500, CV = 1)
  idm.lastnegvisit <- idm(formula01 = Hist(time = list(idm.l.lastnegvisit, idm.r.lastnegvisit), 
                                           event = persist.endpoint) ~ X, 
                          formula02 = Hist(time = obs.dt, event = obs.ds) ~ X,
                          formula12 = Hist(time = obs.dt, event = obs.ds) ~ X,
                          data = as.data.frame(d), maxiter = 500, CV = 1)
  
  # parameter estimates for X
  W.X <- lapply(list(true.01, W.lastvisit, W.last2visit, W.lastnegvisit, W.deathtime), function(x){
    if(!is.null(x$fit)) {x$fit$coef["X",]}else{x$coef["X", ]}
  })
  ic.X <- lapply(list(ic.lastvisit, ic.last2visit, ic.lastnegvisit, ic.deathtime), function(x){
    c(x$coefficients["X"], sqrt(x$var["X","X"]))
  })
  idm.X <- lapply(list(idm.lastvisit, idm.last2visit, idm.lastnegvisit), function(x){
    c(x$coef[1], x$se[1])
  })
  beta01.results <- data.frame(names = c("true.w",
                                         "W.lastvist",
                                         "W.last2visit",
                                         "W.lastnegvisit",
                                         "W.deathtime",
                                         "ic.lastvisit",
                                         "ic.last2visit",
                                         "ic.lastnegvisit",
                                         "ic.deathtime",
                                         "idm.lastvist",
                                         "idm.last2visit",
                                         "idm.lastnegvisit"),
                               do.call(rbind, map(list(W.X, ic.X, idm.X), ~ do.call(rbind, .x))))
  
  # parameter for weib distribution: lambda, gamma
  get.scale <- function(scale, shape){c(shape/(scale^shape), shape)}
  
  W.pars <- lapply(list(W.lastvisit, W.last2visit, W.lastnegvisit, W.deathtime), function(x){
    x$fit$coef[c("lambda", "gamma"), 1]
  })
  ic.pars <- lapply(list(ic.lastvisit, ic.last2visit, ic.lastnegvisit, ic.deathtime), function(x){
    get.scale(exp(x$coefficients["log_scale"]), exp(x$coefficients["log_shape"]))
    
  })
  idm.pars <- lapply(list(idm.lastvisit, idm.last2visit, idm.lastnegvisit), function(x){
    c(x$modelPar[2]^x$modelPar[1], x$modelPar[1])
  })
  
  Wpar.result <- data.frame(names = beta01.results$names[-1],   
                            do.call(rbind, list(W.pars, ic.pars, idm.pars) %>% 
                                      map(~do.call(rbind,.x)) ))
  
  # number of excluded subjects due to 0 survival time
  rn0 <-  rbind(overall = sum(d$tstop == 0),
                last.visit = W.lastvisit$n0,
                last.second.visit = W.last2visit$n0,
                last.neg.test = W.lastnegvisit$n0,
                cens.dead = W.deathtime$n0
                
  )
  rn0 <- data.frame(names = rownames(rn0), n0 = rn0)
  rownames(rn0) <- NULL
  
  # number of events
  rN <- data.frame(n.event = sum(d$ills),
                   n.death = sum(d$ds),
                   n.persist = sum(d$persist.endpoint),
                   n.death.obs = sum(d$obs.ds), # death or censored before the last visit
                   n.event.obs = sum(d$test.result),
                   n.lastnegtest = sum(d$test.result == 1 & d$persist.endpoint == 0),
                   n.lastnegtest.fp = sum(d$test.result == 1 & d$illt > d$tstop)
  )
  
  
  return(list(paraX = beta01.results,
              paraW = Wpar.result,
              n0 = rn0,
              events = rN))
}


###################################### 
# Function to do n simulation run
###################################### 

simulation_nrun <- function(Niter = 100,
                            n,
                            beta01, 
                            beta02, 
                            beta12,
                            shape01,
                            shape02,
                            shape12,
                            scale01,
                            scale02,
                            scale12,
                            cova.prob, 
                            inspectLength, 
                            inspections, 
                            false.positive.rate,
                            min.cens = 0,
                            max.cens = 100,
                            max.int = 500){
  r <- lapply(1:Niter, function(x){
    
    simulation_result <- simulation_1run(n = n,
                                         beta01 = beta01, 
                                         beta02 = beta02, 
                                         beta12 = beta12,
                                         shape01 = shape01, 
                                         shape02 = shape02, 
                                         shape12 = shape12,
                                         scale01 = scale01, 
                                         scale02 = scale02, 
                                         scale12 = scale12,
                                         cova.prob = cova.prob, 
                                         inspectLength = inspectLength, 
                                         inspections = inspections, 
                                         false.positive.rate = false.positive.rate,
                                         min.cens = min.cens, 
                                         max.cens = max.cens, 
                                         max.int = max.int)
    simulation_result
  })
  
  paraX <- do.call(rbind, lapply(r, function(x) x$paraX))
  paraW <- do.call(rbind, lapply(r, function(x) x$paraW))
  n0 <- do.call(rbind, lapply(r, function(x) x$n0))
  nevents <- do.call(rbind, lapply(r, function(x) x$events))
  
  return(list(paraX = paraX,
              paraW = paraW,
              n0 = n0,
              events = nevents))
} 

############################################################## 
# Function to do n simulation run
# output: 
# - Number of events when false positive rate = 1
############################################################## 

fp_nrun <- function(Niter = 100,
                    n,
                    beta01, 
                    beta02, 
                    beta12,
                    shape01,
                    shape02,
                    shape12,
                    scale01,
                    scale02,
                    scale12,
                    cova.prob, 
                    inspectLength, 
                    inspections, 
                    false.positive.rate,
                    min.cens = 0,
                    max.cens = 100,
                    max.int = 500){
  r <- lapply(1:Niter, function(x){
    
    d <- sim.illdeath.Interval.FP(n = n, 
                                  shape01 = shape01, shape02 = shape02, shape12 = shape12,
                                  scale01 = scale01, scale02 = scale02, scale12 = scale12, 
                                  beta01 = beta01, beta02 = beta02, beta12 = beta12, 
                                  cova.prob = cova.prob, inspections = inspections, 
                                  inspectLength = inspectLength,
                                  false.positive.rate =  false.positive.rate, 
                                  min.cens = min.cens, max.cens = max.cens, 
                                  max.int = max.int)
    
    rN <- data.frame(n.event = sum(d$persist.endpoint == 1),
                     n.event.incorrect.time = sum(d$persist.endpoint == 1 & d$illt > d$tstop),
                     n.event.would.not.obverve = sum(d$persist.endpoint == 1 & d$illt > d$tstop & d$n.inspect == d$nrows + 1))
    return(list(events = rN))
  })
  
  
  nevents <- do.call(rbind, lapply(r, function(x) x$events))
  
  return(list(events = nevents))
} 

