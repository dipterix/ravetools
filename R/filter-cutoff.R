# butter_cutoff <- function(type, n, w, r) {
#   # DIPSAUS DEBUG START
#   # w <- 60 / 250
#   # type <- "high"
#   # w <- c(1, 60) / 250
#   # type <- "pass"
#   # r <- 6
#   # n <- 6
#
#   Ws <- sort(w)
#   Rs <- r
#   Rc <- 10 * log10(2)
#
#   # Wpw <- (2/T) * tan(pi * Wp/T)
#   Wsw <- tan(pi * Ws / 2)
#
#   qs <- log(10 ^ (Rs / 10) - 1)
#   qc <- log(10 ^ (Rc / 10) - 1)
#
#   # n <- 0.5 * (qs - qp)/log(ws/wp)
#   # log(ws/wc) =
#   log_ws_wc <- (qs - qc) / ( 2 * n )
#   ws_wc <- exp( log_ws_wc )
#
#   switch(
#     type,
#     "low" = {
#       Wcw <- Wsw / ws_wc
#     },
#     "high" = {
#       # wc <- Wsw
#       # ws <- Wcw
#       Wcw <- Wsw * ws_wc
#     },
#     "pass" = {
#       w02 <- Wsw[[1]] * Wsw[[2]]
#       # wp <- Wpw[2] - Wpw[1]
#       ws <- Wsw[2] - Wsw[1]
#
#       # ws <- ws/wp
#       # wp <- 1
#       wc <- ws / ws_wc
#
#       # solve:
#       # Wcw[2] - Wcw[1] = wc
#       # Wcw[2] * Wcw[1] = w02
#       #
#       # i.e.
#       # Wcw[1]^2 + wc * Wcw[1] - w02 = 0
#       Wcw1 <- (-wc + sqrt(wc^2 + 4 * w02)) / 2
#       Wcw2 <- Wcw1 + wc
#       # 0.005631939 0.441716719
#       Wcw <- c(Wcw1, Wcw2)
#     },
#     "stop" = {
#       w02 <- Wsw[1] * Wsw[2]
#
#       # wp <- w02/(Wpw[2] - Wpw[1])
#       # ws <- w02/(Wsw[2] - Wsw[1])
#       # ws <- ws/wp
#       # wp <- 1
#
#       ws <- w02/(Wsw[2] - Wsw[1])
#       wc <- ws / ws_wc
#       wc <- w02 / wc
#
#       # solve:
#       # Wcw[2] - Wcw[1] = wc
#       # Wcw[2] * Wcw[1] = w02
#       #
#       # i.e.
#       # Wcw[1]^2 + wc * Wcw[1] - w02 = 0
#       Wcw1 <- (-wc + sqrt(wc^2 + 4 * w02)) / 2
#       Wcw2 <- Wcw1 + wc
#       Wcw <- c(Wcw1, Wcw2)
#     }
#   )
#
#   Wc <- atan( Wcw ) * 2 / pi
#
#   # filter <- butter(n, Wc, type)
#   # hw <- freqz2(filter$b, filter$a, fs = 500)
#   # print(hw, cutoffs = -Rs)
#   # filtmag_db(filter$b, filter$a, Wc) + Rc
#   # filtmag_db(filter$b, filter$a, Ws) + Rs
#   Wc
#
# }
#
# cheb1_cutoff <- function(type, n, w, r) {
#   sort(w)
# }
#
# cheb2_cutoff <- function(type, n, w, r) {
#   # DIPSAUS DEBUG START
#   # w <- 60 / 250
#   # type <- "high"
#   # w <- c(1, 60) / 250
#   # type <- "pass"
#   # w <- 0.1
#   # type = "low"
#   # r <- 6
#   # n <- 6
#
#   Ws <- sort(w)
#   Rs <- abs(r)
#   Rc <- 10 * log10(2)
#
#   if( abs(Rs - Rc) < 1e-7 ) {
#     return(Ws)
#   }
#
#   # Wpw <- (2/T) * tan(pi * Wp/T)
#   Wsw <- tan(pi * Ws / 2)
#
#   # Wa <- ws/wp
#   qs <- 10 ^ (abs(Rs) / 10)
#   qc <- 10 ^ (abs(Rc) / 10)
#   # n <- ceiling(acosh(sqrt((cutt_atten - 1)/(pass_atten - 1)))/acosh(Wa))
#   if( qc > qs ) {
#     ws_wc <- cosh( acosh(sqrt((qc - 1)/(qs - 1))) / n )
#   } else {
#     ws_wc <- 1 / cosh( acosh(sqrt((qs - 1)/(qc - 1))) / n )
#   }
#
#
#   switch(
#     type,
#     "low" = {
#       Wcw <- Wsw / ws_wc
#     },
#     "high" = {
#       # wc <- Wsw
#       # ws <- Wcw
#       Wcw <- Wsw * ws_wc
#     },
#     "pass" = {
#       w02 <- Wsw[[1]] * Wsw[[2]]
#       # wp <- Wpw[2] - Wpw[1]
#       ws <- Wsw[2] - Wsw[1]
#
#       # ws <- ws/wp
#       # wp <- 1
#       wc <- ws / ws_wc
#
#       # solve:
#       # Wcw[2] - Wcw[1] = wc
#       # Wcw[2] * Wcw[1] = w02
#       #
#       # i.e.
#       # Wcw[1]^2 + wc * Wcw[1] - w02 = 0
#       Wcw1 <- (-wc + sqrt(wc^2 + 4 * w02)) / 2
#       Wcw2 <- Wcw1 + wc
#       # 0.005631939 0.441716719
#       Wcw <- c(Wcw1, Wcw2)
#     },
#     "stop" = {
#       w02 <- Wsw[1] * Wsw[2]
#
#       # wp <- w02/(Wpw[2] - Wpw[1])
#       # ws <- w02/(Wsw[2] - Wsw[1])
#       # ws <- ws/wp
#       # wp <- 1
#
#       ws <- w02/(Wsw[2] - Wsw[1])
#       wc <- ws / ws_wc
#       wc <- w02 / wc
#
#       # solve:
#       # Wcw[2] - Wcw[1] = wc
#       # Wcw[2] * Wcw[1] = w02
#       #
#       # i.e.
#       # Wcw[1]^2 + wc * Wcw[1] - w02 = 0
#       Wcw1 <- (-wc + sqrt(wc^2 + 4 * w02)) / 2
#       Wcw2 <- Wcw1 + wc
#       Wcw <- c(Wcw1, Wcw2)
#     }
#   )
#
#   Wc <- atan( Wcw ) * 2 / pi
#
#   # filter <- butter(n, Wc, type)
#   # hw <- freqz2(filter$b, filter$a, fs = 500)
#   # print(hw, cutoffs = -Rs)
#   # filtmag_db(filter$b, filter$a, Wc) + Rc
#   # filtmag_db(filter$b, filter$a, Ws) + Rs
#   Wc
#
# }
#
# 0.003660954 0.217171785
# cheb2ord(c(1,50)/250, c(0.5,55)/250, 0.05, 20)
#
# cheb2_cutoff("pass", 12, c(1,50)/250, 20)
#
# cheb2_cutoff("low", 3, 0.1, 6)
# with(
#   gsignal:::cheby2.default(n = 3, Rs = 3, w = cheb2_cutoff("low", 3, 0.1, 6), type = "low"),
#   {
#     check_filter(b,a,cheb2_cutoff("low", 3, 0.1, 6),-3)
#   }
# )
#
# with(
#   gsignal:::cheby2.default(n = 3, Rs = 6, w = 0.1, type = "low"),
#   {
#     check_filter(b,a, cheb2_cutoff("low", 3, 0.1, 6),-3)
#   }
# )
