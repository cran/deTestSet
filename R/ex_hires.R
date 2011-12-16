## =============================================================================
##
## Hires problem, plant physiology
##
## High irradiance response of morphogenesis on the basis of 
## phytochrome; 
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 8
##
## =============================================================================

hires <- function(times = seq(0, 5, by = 0.01), yini = NULL,
                  parms = list(),  ...) {

### derivative function
  hires <- function(t,y,parms) {
    with (as.list(c(y,parms)),{
      dPr     <- -k1*Pr  +  k2*Pfr     + k6*PrX + Oks
      dPfr    <-  k1*Pr  - (k2+k3)*Pfr
      dPrX    <- -(k6+k1)*PrX + k2*PfrX     + k5*PrX2
      dPfrX   <-  k3*Pfr + k1*PrX      -(k4+k2)*PfrX
      dPrX2   <- -(k5+k1)*PrX2 + k2*(PfrX2+PfrX2E)
      dPfrX2  <- -k7*PfrX2*E + k8*PfrX + k1*PrX2 - k2*PfrX2+  k8*PfrX2E
      dPfrX2E <-  k7*PfrX2*E - (k2+k8+k9)*PfrX2E
      dE      <- -k7*PfrX2*E + (k2+k8+k9)*PfrX2E

    list(c(dPr, dPfr, dPrX, dPfrX, dPrX2, dPfrX2, dPfrX2E, dE))
  })
}

### check input 
   parameter <- c(k1 = 1.71, k2 = 0.43, k3 = 8.32, k4 = 0.69, k5 = 0.035,
       k6 = 8.32, k7 = 280, k8 = 0.69, k9 = 0.69, Oks = 0.0007)
   parameter <- overrulepar(parameter, parms, 10)

   if (is.null(yini)) yini <- c(Pr=1, Pfr=0, PrX=0, PfrX=0, 
        PrX2=0, PfrX2=0, PfrX2E=0, E=0.0057)
   checkini(8, yini)


### solve
   out <- ode(func = "hiresfun",  dllname = "deTestSet",
              initfunc = "hirespar", parms = parameter,
              y = yini, times = times, ...)

   return(out)
}

