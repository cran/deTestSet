
### ============================================================================
###      Fekete problem (in stabilized index 2 formulation)
###      index 2 DAE of dimension 160
### ============================================================================

fekete <- function(times = seq(0, 20, by = 0.1), yini = NULL, dyini = NULL,
                   parms = list(), method = "mebdfi", maxsteps = 1e5, ...) {

# No parameters 

### check input 
    if (is.null(yini) | is.null(dyini)) { 
      # Initial conditions in a fortran function
      Init <- .Fortran("fekinit",N=as.integer(160),
                  T=as.double(0),Y=as.double(rep(0.,160)),
                  dY=as.double(rep(0.,160)))
      if (is.null(yini)) yini   <- Init$Y
      if (is.null(dyini)) dyini  <- Init$dY
    }  

    checkini(160, yini, dyini)
    if (is.null(names(yini)))
      names(yini) <- c(
        paste("px",1:20,sep=""),paste("py",1:20,sep=""),paste("pz",1:20,sep=""),
        paste("qx",1:20,sep=""),paste("qy",1:20,sep=""),paste("qz",1:20,sep=""),
        paste("lambda",1:20,sep=""),paste("mu",1:20,sep=""))

### solve
   ind  <- c(6*20,2*20,0)    #  index of the system

   useres <- FALSE
   if (is.character(method)) {
    if (method %in% c("mebdfi", "daspk"))
      useres <- TRUE
   } else  if("res" %in% names(formals(method)))
      useres <- TRUE

    if (useres)
       return( dae(y = yini, dy = dyini, times = times, res = "fekres",
                   nind = ind, method = method,
                   dllname = "deTestSet", initfunc = NULL,
                   parms = NULL, maxsteps = maxsteps))

   fekete <- dae(y = yini, times = times, nind = ind,
          func = "fekfunc", mass = c(rep(1, 120), rep(0, 40)),
          massup = 0, massdown = 0,
          dllname = "deTestSet", initfunc = NULL,
          parms = NULL, method = method,
          maxsteps = maxsteps, ...)

   return(fekete)
}
