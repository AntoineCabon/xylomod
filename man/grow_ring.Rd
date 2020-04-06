\name{grow_ring}
\alias{grow_ring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ring growth simulation
}
\description{
Function to simulate the growth of a ring radial file, including cell number increment by cambial division based on \code{\link{divide}} and increase in cell diameter during cell differentiation based on \code{\link{expand}}.
}
\usage{
grow_ring(ring, data,
          Nc = 8.85, phi0 = 0.13, pi0 = -0.8, CRD0 = 8.3,
          Y_P = 0.05, Y_T = 8, h = 0.043 * 1.8, s = 1.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ring}{
A \code{ring} object, analogous to that returned by the function \code{initialize_ring}
}
  \item{data}{
A dataframe with a column "date" containing \code{Date} class calandar dates, a column "psi" containing water potential in MPa and a column "Tc" containing temperature in °C
}
  \item{Nc}{
The number of cells able to divide in the cambium
}
\item{phi0}{Initial value of cell extensibility (in MPa-1 day-1)}
\item{pi0}{Initial value of cell osmotic potential (in MPa)}
\item{CRD0}{Initial value of cell radial diameter}
\item{Y_P}{Turgor pressure yield threshold (in MPa)}
\item{Y_T}{Temperature yield threshold (in °C)}
\item{h}{Cell wall hardening coefficient (in day-1)}
\item{s}{Cell wall softening coefficient (unitless)}
}
\details{

}
\value{
Returns a \code{ring} object with updated elements:

\item{$divisions}{New lines are appended to the dataframe supplied in the input \code{ring} object for each date supplied in the in the input \code{data}, with the correpsonding daily cell production "P" calculated based on the function \code{divide}.}

\item{$cells}{Add new cells (i.e. lines) and provides the final value of the cell expansion state variables "phi", "pi" and "CRD" for previously present cells and newly added cells, based on the function \code{expand}.

If \code{attr(ring, "cell_wise") == F}, new lines are appended to the dataframe supplied in the input \code{ring} object for each date supplied in the in the input \code{data}. In this case, lines therefore represent theoretical cells formed on each day of the simulation. The contribution of each theoretical cell to the final ring width is thus given by the corresponding daily production "P" and total ring width is obtained as \code{sum(ring$cells$CRD * ring$divisions$P)}.

If \code{attr(ring, "cell_wise") == T}, new lines are appended to the dataframe supplied in the input \code{ring} object for each production of a whole new cell, as calculated based on \code{divide}. In this case, lines therefore cells that were actually simulated to be formed. Total ring width is thus given by \code{sum(ring$cells$CRD)}.
}

\item{cells_historic}{TO BE DONE}

%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cabon A, Fernández‐de‐Uña L, Gea‐Izquierdo G, Meinzer FC, Woodruff DR, Martínez‐Vilalta J, De Cáceres M. 2020a. Water potential control of turgor‐driven tracheid enlargement in Scots pine at its xeric distribution edge. New Phytologist 225: 209–221.

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
\code{\link{initialize_ring}}, \code{\link{expand}}, \code{\link{expand_ring}}, \code{\link{divide}}
}
\examples{
# Example data set with daily temperature and water potential series
data(example_input)

# Select the first year
data <- example_input[1:365,]
data$Tc <- data$Temp_air_mean


### Using cell_wise == T
# Initialize empty ring object with cell wise representation and cell expansion historic
ring <- initialize_ring(NA, cell_wise = T, historic = T)

# Simulate ring growth
ring <- grow_ring(ring, data, Nc = 4)

# Plot the results:
# Simulated cell production rate and total annual cell production
plot(ring$divisions, type = "l", ylab = "Cell production rate (cell day-1)", main = paste0("Annual cell production = ", floor(sum(ring$divisions$P))))

# Cumulative cell production
plot(cumsum(ring$divisions$P)~ring$divisions$date, type = "l", xlab = "date", ylab = "Cell number")

# Simulated cell expansion dynamics
matplot(ring$cells_historic$CRD, type = "l", xlab = "Day of the year", ylab = "Cell radial diameter (micrometers)")

# Resulting tree ring structure
plot(ring$cells$CRD, type = "b", xlab = "Cell index", ylab = "Cell radial diameter (micrometer)")

# Ring width increment
plot(rowSums(ring$cells_historic$CRD, na.rm = T), type = "l", xlab = "Day of the year", ylab = "Ring width (micrometers)")
# Comparison with cell number increment
lines(floor(cumsum(ring$divisions$P))*mean(ring$cells$CRD), col = 2, lty = 2)


### Using cell_wise == F
# Initialize empty ring object with cell wise representation and cell expansion historic
ring <- initialize_ring(NA, cell_wise = F, historic = T)

# Simulate ring growth
ring <- grow_ring(ring, data, Nc = 4)

# Plot the results:
# Simulated cell expansion dynamics
matplot(ring$cells_historic$CRD, type = "l", xlab = "Day of the year", ylab = "Cell radial diameter (micrometers)")

# Resulting tree ring structure
plot(ring$cells$CRD~cumsum(ring$divisions$P), type = "l", xlab = "Cell index", ylab = "Cell radial diameter (micrometer)")

# Ring width increment
plot(rowSums(ring$cells_historic$CRD*matrix(ring$divisions$P, nrow = 365, ncol = 365, byrow = T), na.rm = T), type = "l", xlab = "Day of the year", ylab = "Ring width (micrometers)")
# Comparison with cell number increment
lines(cumsum(ring$divisions$P)*weighted.mean(ring$cells$CRD, ring$divisions$P), col = 2, lty = 2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~grow }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~ring }% __ONLY ONE__ keyword per line
