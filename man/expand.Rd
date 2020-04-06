\name{expand}
\alias{expand}
\alias{expand_ring}
\title{Tracheid expansion model}

\description{
Functions to perform turgor-driven tracheid expansion simulations based on and water potential and temperature inputs.
\code{expand} takes water potential and temperature series inputs and calculates daily expansion dynamics for one cell or a sequence of cell.
\code{expand_ring} takes water potential and temperature series inputs and calculates daily expansion dynamics for cells within a \code{ring} object.
}

\usage{
expand(data,
  phi0 = 0.13, pi0 = -0.8, CRD0 = 8.3,
  Y_P = 0.05, Y_T = 5, h = 0.077, s = 1.8)

expand_ring(ring, data,
  phi0 = 0.13, pi0 = -0.8, CRD0 = 8.3,
  Y_P = 0.05, Y_T = 5, h = 0.077, s = 1.8)
}

\arguments{
\item{ring}{A \code{ring} object as returned by the function \code{initialize_ring}}
\item{data}{A dataframe containing a column "date" with dates specified using \code{Date} class, a column "psi" containing water potentials in MPa and a column "Tc" containing temperature in °C}
\item{phi0}{Initial value of cell extensibility (in MPa-1 day-1)}
\item{pi0}{Initial value of cell osmotic potential (in MPa)}
\item{CRD0}{Initial value of cell radial diameter}
\item{Y_P}{Turgor pressure yield threshold (in MPa)}
\item{Y_T}{Temperature yield threshold (in °C)}
\item{h}{Cell wall hardening coefficient (in day-1)}
\item{s}{Cell wall softening coefficient (unitless)}
}

\value{
\code{expand} returns a dataframe containing the cell state variables ("phi", "pi" and "CRD", in columns) for each date of input data (in rows).

\code{expand_ring} returns a ring object containing an updated "cells" dataframe with the value of the cell's state variable at the end of each simulation. If \code{attr(ring, "historic") == TRUE}, then the dynamics of the state variables during the simulation are also returned in the list \code{cells_historic} containing one matrix for each state variable with cells in column and dates in rows.
}
\seealso{
\code{\link{divide}}, \code{\link{grow_ring}}, \code{\link{initialize_ring}}
}

\examples{
# Example data set containing air temperature and soil water potential daily measurements over five years
data("example_input")

# Select 100 days starting at the beginning of the growing season (DOY 150)
data <- example_input[150:250,]
data$Tc <- data$Temp_air_mean

### Calculate cell expansion for a cell formed on DOY 150
expansion <- expand(example_input)

# Plot the result
plot(CRD~date, data = expansion, type = "l", y
     lab="Cell radial diameter (micrometers)")

### Calculate cell expansion for three different cells formed at specified dates
# Initialize ring object
formation_dates <- c(as.Date("2010-06-01"),
                     as.Date("2010-06-15"),
                     as.Date("2010-06-30"))
ring <- initialize_ring(formation_dates, historic = T)

# Calculate cell expansion
expand_ring(ring, data)

# Print cells state at the end of simulation
ring$cells

# Plot the historic of cell radial diameter variations during the simulation
matplot(ring$cells_historic$CRD, xlab = "Time step",
        ylab = "Cell radial diameter (micrometer)")

}

\references{
Cabon A, Fernández‐de‐Uña L, Gea‐Izquierdo G, Meinzer FC, Woodruff DR, Martínez‐Vilalta J, De Cáceres M. 2020a. Water potential control of turgor‐driven tracheid enlargement in Scots pine at its xeric distribution edge. New Phytologist 225: 209–221.

}
\author{
Antoine Cabon, CTFC
}
