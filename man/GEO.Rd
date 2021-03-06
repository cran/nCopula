% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Structure_Distribution_GEO.R
\name{GEO}
\alias{GEO}
\title{Construction of a GEO Mother or Child Class Object}
\usage{
GEO(par, unif, structure)
}
\arguments{
\item{par}{parameter of the distribution.}

\item{unif}{uniform structure, a numeric vector of grouped.
numbers, i.e. c(1,2,3) is translated as being c(u1, u2, u3).}

\item{structure}{nesting structure of the form

X(par1, c(i,...), list(Y(par2, c(j,...), NULL),
                       Z(par3, c(k,...), NULL))),

where X, Y, and Z are compatible functions (see 'details').
It is to note that if structure is NULL, the function will automatically
be of class Child. For continuous distributions (i.e. GAMMA), structure is
always NULL.}
}
\description{
Constructs either a GEO Mother or Child class object for
a given parameter, arguments, and nesting structure.
}
\examples{
GEO(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
                    GEO(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
                                        GAMMA(1/30, c(3,4), NULL)))))
}
\seealso{
Other mother or child class objects.: \code{\link{GAMMA}},
  \code{\link{LOG}}
}
\author{
Simon-Pierre Gadoury
}
