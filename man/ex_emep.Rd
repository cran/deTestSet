\name{emep}
\alias{emep}
\title{Emep MSC-W Ozone Chemistry Problem, ODE}
\description{
The problem is a stiff system of 66 ordinary differential equations. 
The 'Mathematics and the Environment'
project group at CWI contributed this problem to the test set. 
The software part of the
problem is in the file emep.f available at [MM08].
}
\usage{
emep (times = seq(14400, 417600, by = 400), yini = NULL, 
      parms = list(), printmescd = TRUE, method = bimd, 
      atol = 0.1, rtol = 1e-5, maxsteps = 1e5, ...)
}

\arguments{
  \item{yini }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time.
  }
  \item{parms }{list of parameters that overrule the default parameter values
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
   \item{method }{the solver to use
  }
   \item{printmescd }{if TRUE  the mixed error significant digits computed using
    the  reference solution at time   417600 are printed
    }
  \item{maxsteps }{maximal number of steps per output interval taken by 
    the solver
  }  
  \item{... }{additional arguments passed to the solver .
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{yini}, plus an additional column (the first) 
  for the time value.

  There will be one row for each element in \code{times} unless the
  solver returns with an unrecoverable error.  If
  \code{yini} has a names attribute, it will be used to label the columns
  of the output value.
      
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>

Francesca Mazzia <mazzia@dm.uniba.it>
}
\examples{
out <- emep()
plot(out, lwd = 2, col = "darkblue", 
  which = c("NO", "NO2", "SO2", "CH4", "O3", "N2O5"))

plot(out, col = "darkblue", lwd = 2, which = 1:16)
mtext(side = 3,line = -1.5, "emep", cex = 1.25, outer = TRUE)

# compare with reference solution (component 36 and 38 not included)       
refsol <- reference("emep")
inderr <- c(1:35,37,39:66)
max(abs(out[nrow(out),inderr+1] - refsol[inderr])/refsol[inderr])     

}
\references{
\url{http://archimede.dm.uniba.it/~testset}

[MM08] F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers, 
release 2.4. Department
of Mathematics, University of Bari and INdAM, Research Unit of Bari, February
2008.  

[SASJ93] D. Simpson, Y. Andersson-Skold, and M.E. Jenkin. Updating the chemical scheme for the
EMEP MSC-W model: Current status. Report EMEP MSC-W Note 2/93, The Norwegian
Meteorological Institute, Oslo, 1993.

[Sim93] D. Simpson. Photochemical model calculations over Europe for two extended summer
periods: 1985 and 1989. model results and comparisons with observations. Atmospheric
Environment, 27A:921-943, 1993.

[VS94] J.G. Verwer and D. Simpson. Explicit methods for stiff ODEs from atmospheric chemistry.
Report NM-R9409, CWI, Amsterdam, 1994.
}
\details{

  The default parameters are:

  c = 1.6e-8 , cs = 2e-12 , cp = 1e-8  , r = 25e3 , rp = 50,
          lh = 4.45   , ls1 = 2e-3 , ls2 = 5e-4 , ls3 = 5e-4,
          rg1 = 36.3  , rg2 = 17.3 , rg3 = 17.3 , ri = 50 , rc = 600,
          gamma = 40.67286402e-9 , delta = 17.7493332
          
  
}
\keyword{utilities}

