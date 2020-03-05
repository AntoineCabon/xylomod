\name{expansion}
\alias{expansion}
\alias{expansion_seq}
\title{Tracheid expansion model}
\usage{
expansion(psi, Tc, start = 1,
  phi0 = 0.13, pi0 = -0.8, CRD0 = 8.3,
  Y_P = 0.05, Y_T = 5, h = 0.077, s = 1.8)
expansion_seq(psi, Tc, start_vec,
  phi0 = 0.13, pi0 = -0.8, CRD0 = 8.3,
  Y_P = 0.05, Y_T = 5, h = 0.077, s = 1.8)
}

\arguments{
\item{psi}{Water potential (in MPa)}
\item{Tc}{Temperature (in deg C)}
\item{start}{An integer scalar giving the vectorial position of cell expansion onset}
\item{start_vec}{An integer vector giving the vectorial position of cell expansion onset}
\item{phi0}{Initial value of cell extensibility (in MPa-1 day-1)}
\item{pi0}{Initial value of cell osmotic potential (in MPa)}
\item{CRD0}{Initial value of cell radial diameter}
\item{Y_P}{Turgor pressure yield threshold (in MPa)}
\item{Y_T}{Temperature yield threshold (in deg C)}
\item{h}{Cell wall hardening coefficient (in day-1)}
\item{s}{Cell wall softening coefficient (unitless)}
}

\description{
Functions to perform turgor-driven tracheid expansion simulations based on and water potential and temperature inputs.
\code{expansion} takes water potential and temperature series inputs and returns daily expansion dynamics for one cell or a sequence of cell.
\code{expansion_seq} takes water potential and temperature series inputs and returns daily expansion dynamics for as many cells as designated by \code{start_vec}.
}
\examples{
load("~/These/Xylogenesis - cell enlargement/RData/met_swc_psi2.RData")
subs <- 1:365
dates <- met$date[subs]
expansion_res <- expansion_seq(psi = met$psi[subs],
  Tc = met$MeanTemperature[subs], start = seq_along(subs), Y_T = 8)

plot(expansion_res[[1]]$V~dates, ylim = c(8,50),type = "l")
for(i in seq(1,length(subs), by = 10)){
 lines(expansion_res[[i]]$V~dates)
}
}
\references{
Cabon A, Fernández‐de‐Uña L, Gea‐Izquierdo G, Meinzer FC, Woodruff DR, Martínez‐Vilalta J, De Cáceres M. 2020. Water potential control of turgor‐driven tracheid enlargement in Scots pine at its xeric distribution edge. New Phytologist 225: 209–221.
}
\author{
Antoine Cabon, CTFC
}
