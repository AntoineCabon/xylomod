#
# load("~/These/Xylogenesis - cell enlargement/RData/met_swc_psi2.RData")
# subs <- which(met$date %in% seq.Date(as.Date("2013-01-01"), as.Date("2013-12-31"), "day"))
# dates <- met$date[subs]
# expansion_res <- expansion_seq(psi = met$psi[subs], Tc = met$MeanTemperature[subs], start = seq_along(subs), Y_T = 8)
#
# # plot(expansion_res[[1]]$CRD~dates, ylim = c(8,50),type = "l")
# # for(i in seq(1,length(subs), by = 10)){
# #   lines(expansion_res[[i]]$CRD~dates)
# # }
#
# division_res <- division(psi = met$psi[subs], Tc = met$MeanTemperature[subs], Y_T = 8, Nc = 1)
# plot(cumsum(division_res$P)~dates, type = "l")
#
# starts <- sapply(expansion_res, FUN = function(x)unique(x$start))
# active_mat <- matrix(seq_along(subs), nrow = length(subs), ncol = length(starts), byrow = F) >= matrix(starts, nrow = length(subs), ncol = length(starts), byrow = T)
# division_mat <- matrix(division_res$P, nrow = length(subs), ncol = length(starts), byrow = T) * active_mat
# expansion_mat <- array(unlist(expansion_res),
#                          dim = c(length(subs), ncol(expansion_res[[1]]), length(subs)),
#                          dimnames = list(Date = subs, Var = colnames(expansion_res[[1]]), Start = starts))[,"CRD",]
# growth_mat <- expansion_mat*division_mat
# growth_vec <- rowSums(growth_mat)
#
# ring <- data.frame(Date = dates, P = division_res$P, N = cumsum(division_res$P),
#                    RW = growth_vec, CRD = expansion_mat[nrow(expansion_mat),])
#
# library(ggplot2)
# ggplot(ring[ring$N>=1,], aes(x = RW, y = CRD))+
#   geom_line()
