\name{divide}
\alias{divide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cambial cell division model
}
\description{
Function to perform simulations of daily cambial cell divisions and wood cell radial file production based on the effect of temperature and and water potential inputs on turgor-driven cell expansion.
}
\usage{
divide(data, Nc = 8.85, phi = 0.13, pi = -0.8, Y_P = 0.05, Y_T = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A dataframe containing a column "date" with dates specified using \code{Date} class, a column "psi" containing water potentials in MPa and a column "Tc" containing temperature in °C
}
  \item{Nc}{
A numeric scalar giving the number of cells able to divide in the cambium (does not need to be integer)
}
  \item{phi}{
Cell wall extensibility (in day-1 MPa-1)
}
  \item{pi}{
Osmotic potential (in MPa)
}
  \item{Y_P}{
Turgor pressure yield threshold (in MPa)
}
  \item{Y_T}{
Temperature yield treshold (in °C)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe containing P, the number of new cells produced by a cambium radial file on each date.
}
\references{
Cabon A, Peters RL, Fonti P, Martínez‐Vilalta J, De Cáceres M. 2020b. Temperature and water potential co‐limit stem cambial activity along a steep elevational gradient. New Phytologist: nph.16456.
}
\author{
Antoine Cabon, CTFC
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{expand}}, \code{\link{grow_ring}, \code{\link{Date}}}
}
\examples{
data("example_input") # example data set containing air temperature and soil water potential daily measurements over five years

data <- example_input[1:365,] # select only the first year
data$Tc <- data$Temp_air_mean

divisions <- divide(data)
plot(P~as.Date(date), data = divisions, type = "l", xlab="Date", ylab="Cell production (cell day-1)")
sum(divisions$P) # gives the total number of cells formed during the year
}
