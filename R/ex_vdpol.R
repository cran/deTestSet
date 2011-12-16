## =============================================================================
##
## van der Pol problem
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 2
##
## =============================================================================

vdpol <- function(times = 0:2000, yini = NULL, parms = list(), ...) {

### derivative function
  Vdpol <- function(t,y,mu) {
    list(c(
    y[2],
    mu * (1 - y[1]^2) * y[2] - y[1]
  ))
}

### check input 
    parameter <- c(mu = 1000)

    parameter <- overrulepar(parameter, parms, 1)

    if (is.null(yini))  yini <- c(y1 = 2, y2 = 0) 
    checkini(2, yini)


### solve
    out <- ode(func = Vdpol, parms = parameter, y = yini, times = times,
       ...)

    return(out)
}
