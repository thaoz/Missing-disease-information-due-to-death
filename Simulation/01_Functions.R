# Auxiliary functions for Death-bias project
# author: Thao Le
# date 13 Dec 2022
# modified date: 17 Feb 2023

################################################################################
# simulate illness death data (clock forward) for weibull distribution 
# (shape = 1 ~ exponential distribution)
################################################################################
# simulate clock forward illness death data
sim.illdeath.fw <- function(n,
                            shape01 = NULL,
                            shape02 = NULL,
                            shape12 = NULL,
                            scale01 = NULL,
                            scale02 = NULL,
                            scale12 = NULL,
                            beta01,
                            beta02,
                            beta12,
                            cova.prob,
                            rndDigits = 3,
                            min.cens = 0,
                            max.cens = 100, 
                            max.int = 500){
  
  # covariate
  cova <- rbinom(n, size = 1, prob = cova.prob)
  
  # survival time 0->1 and 0->2
  stime = numeric(n)
  inc = 1
  while (inc <= n) {
    u <- runif(1)
    cova_i = cova[[inc]]
    temp <- function(x,y){
      scale01*exp(beta01*cova_i)*x^shape01 + scale02*exp(beta02*cova_i)*x^shape02 + y
    }
    if (temp(0, log(1 - u)) * temp(max.int, log(1 - u)) < 0) {
      res <- uniroot(temp, c(0, max.int), tol = 0.0001, y = log(1 - u))
      if(res$root>0){
        stime[inc] <- res$root
        inc <- inc + 1
      }else{
        inc <- inc
      }
      
    }
    else {
      inc <- inc
    }
    
  }
  
  # cause-specfic hazards
  h01 = shape01*scale01*stime^(shape01-1)*exp(beta01*cova)
  h02 = shape02*scale02*stime^(shape02-1)*exp(beta02*cova)
  h0 = h01 + h02
  
  # cause
  f.cause = rbinom(n, size = 1, prob = h01/h0) # cause 1, 0
  f.cause <- ifelse(f.cause == 1, 1, 2) # cause 1, 2
  
  ## Transition 12: Clock forward ##
  ## Generate time from a conditional distribution ## 
  
  # generate waiting time for those from 1
  draw.12 <- function(n, shape12, scale12, T1, cova1){
    p <- runif(n)
    
    ( T1^shape12 - log(1-p)/(scale12*exp(beta12*cova1)))^(1/shape12) 
    
  }
  n12 <- sum(f.cause == 1)
  cova12 <- cova[f.cause == 1]
  
  t12 = draw.12(n = n12,
                shape12 = shape12, 
                scale12 = scale12, 
                T1 = stime[f.cause==1], 
                cova1 = cova12)
  sstime = stime
  sstime[f.cause == 1] <- t12
  illt = stime
  ills = as.numeric(f.cause == 1)
  dt = sstime
  ds = 1
  
  # create right censored time
  cens.times <- runif(n, min.cens, max.cens)
  
  # death time and status, after right censoring
  ds = ds*as.numeric(dt<=cens.times)
  dt = pmin(dt, cens.times)
  
  ills = ills*as.numeric(illt<=cens.times)
  illt = pmin(illt, cens.times)
  s.data = data.frame(id =1:n,
                      X = cova,
                      illt = illt,
                      dt = dt,
                      ills = ills,
                      ds = ds)
  return(s.data)
}
################################################################################
# simulate illness death data with interval censored event time for 
# weibull distribution (shape = 1 ~ exponential distribution)
################################################################################

sim.illdeath.Interval.FP <- function(n,
                                     shape01 = NULL,
                                     shape02 = NULL,
                                     shape12 = NULL,
                                     scale01 = NULL,
                                     scale02 = NULL,
                                     scale12 = NULL,
                                     beta01,
                                     beta02,
                                     beta12,
                                     cova.prob,
                                     inspections,
                                     inspectLength,
                                     false.positive.rate = 0.1,
                                     min.cens = 0,
                                     max.cens = 100, 
                                     max.int = 500){
  
  # simulate complete illness death data
  d <- sim.illdeath.fw(n = n,
                       shape01 = shape01,
                       shape02 = shape02,
                       shape12 = shape12,
                       scale01 = scale01,
                       scale02 = scale02,
                       scale12 = scale12,
                       beta01 = beta01,
                       beta02 = beta02,
                       beta12 = beta12,
                       cova.prob = cova.prob,
                       min.cens = min.cens,
                       max.cens = max.cens, 
                       max.int = max.int)
  # simulate false positive time
  if (false.positive.rate >0){
    fptime.d <- data.frame(id = 1:n,
                           fptime = rexp(n = n, rate = false.positive.rate))
    
  }else{
    fptime.d <- data.frame(id = d$id,
                           fptime = d$illt)
  }
  
  d <- merge(d, fptime.d, by = "id")
  
  # simulate interval censored data
  
  get.ic.time.fp <- function(illt, ills, dt, id, fptime){
    tstop = test.result = 0
    i = 1
    tstart = 0
    .fptime = fptime
    while(tstart < dt & i <= inspections ){
      tstop[i] = tstart + inspectLength + runif(1, 0, 0.01)
      if(tstop[i] < dt){
        if(.fptime >= illt){
          test.result[i] = as.numeric(tstop[i]>=illt)*ills
        }else{
          test.result[i] = as.numeric(tstop[i] >= .fptime)
          .fptime = ifelse(test.result[i] == 1, illt, .fptime)
        }
        tstart = tstop[i]
        i = i + 1
      }else{
        tstop[i] = tstart
        break
      }
      
      
    }
    
    tstart = c(0, tstop)
    n.inspect = i - 1
    max.length = length(test.result)
    max.fu = max(inspections*(inspectLength + runif(1,0,0.01)), max(tstop))
    
    tmp <- data.frame(id = id,
                      n.inspect = n.inspect,
                      ills = ills,
                      illt = illt,
                      fptime = fptime,
                      dt = dt,
                      ds = 1,
                      max.fu = max.fu,
                      max.length = max.length,
                      tstart = tstart[1:max.length],
                      tstop = tstop[1:max.length],
                      test.result = test.result)
    return(tmp)
    
  }
  d.ic <- pmap_df(.l = list(dt = d$dt, 
                                     illt = d$illt, 
                                     ills = d$ills,
                                     id = d$id,
                                     fptime = d$fptime), 
                           .f = get.ic.time.fp)
  # max.fu = inspections*(inspectLength + runif(1,0,0.01))
  # derive persistent endpoint and obtain pat level data
  d.ic.pat <- d.ic %>% 
    group_by(id) %>% 
    mutate(persist = as.numeric(test.result == 1 & lag(test.result, default = 0) == 1),
           persist.endpoint = ifelse(lead(persist, default = 0) == 1,1,0),
           included.row = (lead(persist, default = 0) == 1| persist == 0 & row_number() == n())) %>% 
    filter(included.row == TRUE) %>% 
    group_by(id) %>% 
    filter(row_number() == 1)
  
  d.ic.fp <- d.ic.pat %>% 
    left_join(d %>% select(id, X), by = "id") %>% 
    mutate(obs.ds = ds*as.numeric(dt <= max.fu),
           obs.dt = pmin(dt, max.fu)) %>% 
    select(id, X, illt, ills, dt, ds, fptime, n.inspect, max.fu, tstart, tstop, test.result, persist.endpoint, obs.ds, obs.dt)
  return(d.ic.fp)
}

