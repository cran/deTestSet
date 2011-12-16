## =============================================================================
##
## nand gate problem,  
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     index 0 IDE of dimension 14       
##
## =============================================================================
      
#-----------------------------------------------------------------------
#
# The network equation describing the nand gate
#             C[Y] * Y' - f[Y,t] = 0
# 
# ---------------------------------------------------------------------
nand <- function(times = 0:80, yini = NULL, dyini = NULL,
                 parms=list(), method = "mebdfi", maxsteps = 1e5, ...) {

### check input 
    parameter <- c(RGS = 4, RGD = 4, RBS = 10, RBD = 10,
       CGS = 0.6e-4, CGD = 0.6e-4, CBD = 2.4e-5, CBS = 2.4e-5,
       C9 = 0.5e-4, DELTA = 0.2e-1, CURIS = 1.e-14, VTH = 25.85,
       VDD = 5., VBB = -2.5 )

    parameter <- overrulepar(parameter, parms, 14)

    if (is.null(yini))
      yini <- with (as.list(parameter),
       c(5,5,VBB,VBB,5,3.62385,5,VBB,VBB,3.62385,0,3.62385,VBB,VBB))
      
    if (is.null(dyini))
      dyini <- rep(0, 14)

    checkini(14, yini, dyini)
    
### solve

   ind <- c(14, 0, 0)  # index of the system

   out <- dae(y = yini, dy=dyini, times = times, res = "nandres", nind = ind,
          dllname = "deTestSet",  initfunc = "nandpar", method = method,
          parms = parameter, maxsteps = maxsteps, ...)

  return(out)
}

