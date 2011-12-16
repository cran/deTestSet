\name{deTestSet-package}
\alias{deTestSet-package}
\alias{deTestSet}
\docType{package}
\title{
  Solvers and Test Set for Initial Value Problems of Ordinary Differential 
  Equations (ODE), Partial Differential Equations (PDE) and for 
  Differential Algebraic Equations (DAE)
}
\description{
  R package deTestSet contains the R-version of the ODE and DAE initial value 
  problems test set from \url{http://www.dm.uniba.it/\~testset}.  

  If the test model problem is small enough, then it is implemented in
  pure R. For larger models, the problem is specified in FORTRAN code.

  These implementations were compiled as DLLs, and included in the package.
  The code of these models can be found in the packages
  \url{inst/doc/examples/dynload} subdirectory.
  
  In addition to all solvers present in package \code{deSolve}, 
  \code{deTestSet} contains the initial value problem solvers 
  \code{gamd}, and \code{mebdfi}, implementing a generalised adams method
  and a differential algebraic equation solver of index up to three.
  

}
\details{
  \tabular{ll}{
    Package: \tab deTestSet\cr
    Type:    \tab Package\cr
    Version: \tab 1.0\cr
    License: \tab GNU Public License 2 or above\cr
  }

 
}
\author{
  Karline Soetaert (Maintainer),
  
  Jeff Cash,
  
  Francesca Mazzia
}
\seealso{
  \code{\link[deSolve]{ode}} for a general interface to most of the ODE solvers
  from package \code{deSolve}
  
  \code{\link[deSolve]{ode.1D}}, \code{\link[deSolve]{ode.2D}}, \code{\link[deSolve]{ode.3D}},
  for integrating 1-D, 2-D and 3-D models from package \code{deSolve}
  
  \code{\link{dae}}, a general interface to the dae solvers, including
  \code{mebdfi}, \code{gamd}, and \code{daspk} and \code{radau} (deSolve)

}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(caraxis)
example(nand)
example(andrews)

## open the directory with R sourcecode examples
browseURL(paste(system.file(package = "deTestSet"), "/doc/examples", sep = ""))
## open the directory with C and FORTRAN sourcecode examples
browseURL(paste(system.file(package = "deTestSet"), "/doc/examples/dynload", sep = ""))

## show package vignette with how to use the test set
## + source code of the vignette
vignette("deTestSet")
edit(vignette("deTestSet"))
}
}

\keyword{ package }