\name{division}
\alias{division}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cambial cell division model
}
\description{
Function to perform simulations of daily cambial cell division and tracheid production based on the effect of temperature and and water potential inputs on turgor-driven cell expansion.
}
\usage{
division(psi, Tc, Nc = 8.85, phi = 0.13, pi = -0.8, Y_P = 0.05, Y_T = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{
Water potential vector (in MPa)
}
  \item{Tc}{
Temperature vector (in dec C)
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
Temperature yield treshold (in deg C)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cabon A, Peters RL, Fonti P, Martínez‐Vilalta J, De Cáceres M. 2020. Temperature and water potential co‐limit stem cambial activity along a steep elevational gradient. New Phytologist: nph.16456.
}
\author{
Antoine Cabon, CTFC
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load("~/These/Xylogenesis - cell enlargement/RData/met_swc_psi2.RData")
subs <- 1:365
dates <- met$date[subs]

division_res <- division(psi = met$psi[subs], Tc = met$MeanTemperature[subs], Y_T = 8, Nc = 1)
plot(cumsum(division_res$P)~dates, type = "l")

## Can be coupled with the expansion model in order to estimate ring growth and structure
expansion_res <- expansion_seq(psi = met$psi[subs], Tc = met$MeanTemperature[subs], start = seq_along(subs), Y_T = 8)

starts <- sapply(expansion_res, FUN = function(x)unique(x$start))
active_mat <- matrix(seq_along(subs), nrow = length(subs), ncol = length(starts), byrow = F) >= matrix(starts, nrow = length(subs), ncol = length(starts), byrow = T)
division_mat <- matrix(division_res$P, nrow = length(subs), ncol = length(starts), byrow = T) * active_mat
expansion_mat <- array(unlist(expansion_res),
                         dim = c(length(subs), ncol(expansion_res[[1]]), length(subs)),
                         dimnames = list(Date = subs, Var = colnames(expansion_res[[1]]), Start = starts))[,"CRD",]
growth_mat <- expansion_mat*division_mat
growth_vec <- rowSums(growth_mat)

ring <- data.frame(Date = dates, P = division_res$P, N = cumsum(division_res$P),
                   RW = growth_vec, CRD = expansion_mat[nrow(expansion_mat),])

library(ggplot2)
ggplot(ring[ring$N>=1,], aes(x = RW, y = CRD))+
  geom_line()
