\name{reference}
\alias{reference}
\title{Reference Value of Test Set Problems}
\description{
Estimates the reference solution of the problem
}
\usage{
reference (name = c("andrews", "beam", "caraxis", "crank", "E5", 
    "emep", "fekete", "vdpol", "hires", "nand", "orego", 
    "pleiades", "pollution", "ring", "rober", "transistor",
    "tube", "twobit", "wheelset")) 
}

\arguments{
  \item{name }{the name of the problem whose reference solution is to be
    estimated
  }
}
\value{
  A vector with the reference solution
      
}
\author{Karline Soetaert 

Francesca Mazzia
}
\examples{
reference("ring")
}
\references{
\url{http://www.dm.uniba.it/~testset}
}
\keyword{utilities}

