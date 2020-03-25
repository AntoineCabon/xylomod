# #
# load("~/Documents/These/Xylogenesis - cell enlargement/RData/met_swc_psi2.RData")
# subs <- which(met$date %in% seq.Date(as.Date("2013-01-01"), as.Date("2013-12-31"), "day"))
# dates <- met$date[subs]
# expand_res <- expand_seq(psi = met$psi[subs], Tc = met$MeanTemperature[subs], start = seq_along(subs), Y_T = 8)
#
# # plot(expand_res[[1]]$CRD~dates, ylim = c(8,50),type = "l")
# # for(i in seq(1,length(subs), by = 10)){
# #   lines(expand_res[[i]]$CRD~dates)
# # }
#
# divide_res <- divide(psi = met$psi[subs], Tc = met$MeanTemperature[subs], Y_T = 8, Nc = 1)
# plot(cumsum(divide_res$P)~dates, type = "l")
#
# starts <- sapply(expand_res, FUN = function(x)unique(x$start))
# active_mat <- matrix(seq_along(subs), nrow = length(subs), ncol = length(starts), byrow = F) >= matrix(starts, nrow = length(subs), ncol = length(starts), byrow = T)
# divide_mat <- matrix(divide_res$P, nrow = length(subs), ncol = length(starts), byrow = T) * active_mat
# expand_mat <- array(unlist(expand_res),
#                          dim = c(length(subs), ncol(expand_res[[1]]), length(subs)),
#                          dimnames = list(Date = subs, Var = colnames(expand_res[[1]]), Start = starts))[,"CRD",]
# growth_mat <- expand_mat*divide_mat
# growth_vec <- rowSums(growth_mat)
# #
# ring <- data.frame(Date = dates, P = divide_res$P, N = cumsum(divide_res$P),
#                    RW = growth_vec, CRD = expand_mat[nrow(expand_mat),])
#
# library(ggplot2)
# ggplot(ring[ring$N>=1,], aes(x = RW, y = CRD))+
#   geom_line()

#######################################
# Test ring_expand

# Initiate ring object
# ring <- list("phi"=matrix(1e-2), "pi"=matrix(-1), "CRD"=matrix(1), "P"=c(0))
#
# expand_ring(ring, -0.5, 20)
# grow_ring(ring, -0.5, 20)
