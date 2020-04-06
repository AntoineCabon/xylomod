\name{initialize_ring}
\alias{initialize_ring}
\alias{ring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialization of \code{ring} object
}
\description{
\code{initialize_ring} configures a \code{ring} readable by \code{_ring} functions from package \code{xylomod}. A \code{ring} object is not a formal class but rather an object structure analogous to that derived by \code{initialize_ring}.
}
\usage{
initialize_ring(formation_date = NA, cell_wise = FALSE, historic = FALSE,
                phi0 = 0.13, pi0 = -0.8, CRD0 = 8.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formation_date}{
A vector containing \code{\link{Date}}s indicating dates of cell formation. If NA is suplied (default) the \code{ring} object is initialize with 0 cells. Else, the number of cells initialized is equal to length of the supplied vector.
}
\item{cell_wise}{
A boolean specifying the behaviour of the \code{cells} slot when using the function \code{grow_ring}. If \code{FALSE} (default), a new "cell", i.e. a line, is added at each time step. Else, cells are added for each whole division, as simulated by \code{divide}.
}
\item{historic}{
A boolean indicating if the historic of cell expansion is to be saved in a \code{cells_historic} slot. if \code{FALSE} (default), only the last state of the cells is returned in the \code{cells} slot.
}
\item{phi0}{
Initial value of cell extensibility (in MPa-1 day-1)
}
\item{pi0}{
Initial value of cell osmotic potential (in MPa)
}
\item{CRD0}{
Initial value of cell radial diameter
}
}

\details{
}

\value{
\code{initialize_ring} return a \code{ring} object, that is a list containing at least two data frames \code{divisions} and \code{cells} and two boolean attributes \code{cell_wise} and \code{historic}, which value is inherited from correponding arguments in \code{initialize_ring}. If \code{historic == TRUE}, then an empty list, "cells_historic" is added to the \code{ring} object. The \code{divisions} dataframe contains a column "date" and a column "P", and the \code{cells} dataframe contains the columns "formation_date", "phi", "pi" and "CRD", with as many rows as the number of dates supplied in the argument \code{formation_date}.
}

\references{
%% ~put references to the literature/web site here ~
}

\author{
Antoine Cabon, CTFC
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{expand_ring}}, \code{\link{grow_ring}}
}
\examples{
initialize_ring() # returns an empty ring object

formation_dates <- c(as.Date("2010-06-01"), as.Date("2010-06-15"), as.Date("2010-06-30"))
initialize_ring(formation_dates) # returns a ring object with three pre-initialized cells with the specified formation dates
}
