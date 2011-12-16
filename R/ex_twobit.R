### ============================================================================
###      two bit adding unit problem 
###      index 1 DAE of dimension 350
### ============================================================================

twobit <- function(times = seq(0, 320, by = 0.5),
   yini = NULL, dyini = NULL,
   method = "radau", atol = 1e-4, rtol = 1e-4, 
   maxsteps = 1e5, hmax = 0.1, ...) {

# No parameters 

### check input 
    N <- 350
    if (is.null(yini) | is.null(dyini)) { 
      # Initial conditions in a fortran function
      Init <- .Fortran("twobinit",N=as.integer(N),
                  T=as.double(0),Y=as.double(rep(0.,N)),
                  dY=as.double(rep(0.,N)))
      if (is.null(yini)) yini   <- Init$Y
      if (is.null(dyini)) dyini  <- Init$dY
    }  

    if (is.null(names(yini))) names(yini) <- 
     c(paste("y",1:175,sep=""),  paste("x",1:175,sep=""))
    checkini(350, yini, dyini)

### solve
   ind  <- c(N,0,0)    #  index of the system
   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))

    if (useres)
      return( dae(y = yini, dy = dyini, times = times,
          res = "twobres", nind = ind, dllname = "deTestSet",
          initfunc = NULL, parms = NULL, method = method,
          atol = atol, rtol = rtol, maxsteps = maxsteps, hmax = hmax, ...))
    
   twobit <- dae(y = yini, dy = dyini, times = times, nind = ind,
          func = "twobfunc", mass =  c(rep(1, 175), rep(0, 350-175)),
          massup = 0, massdown = 0,
          dllname = "deTestSet", initfunc = NULL,
          parms = NULL, method = method,
          atol = atol, rtol = rtol, maxsteps = maxsteps, hmax = hmax, ...)

   return(twobit)
}
