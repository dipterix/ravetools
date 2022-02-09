require(ravetools)
x <- rnorm(2000*60*5)

list2env(list(freqs = seq(2,200,4), srate = 2000, wave_num = c(2,20), trend = "constant"), .GlobalEnv)
ratio = (log(max(wave_num)) - log(min(wave_num))) / (log(max(freqs)) - log(min(freqs)))
wave_num = round(exp((log(freqs) - log(min(freqs))) * ratio + log(min(wave_num))))

gc()
profile <- function(expr, env = parent.frame()){
  expr <- substitute(expr)
  gc()
  m1 <- lobstr::mem_used()
  re <- system.time({
    eval(expr, envir = env)
  }, gcFirst = FALSE)
  m2 <- lobstr::mem_used()
  gettextf("Total elapsed time: %.2f, memory changed: %.2f MB", re["elapsed"], (m2 - m1) / 1024^2)
}
# profile <- profvis::profvis

# gctorture2(1000, inhibit_release = TRUE); gc()
profile({
  y2 <- morlet_wavelet(x, freqs, srate, wave_num, trend = trend, precision = "float")
  gc()
})

profile({
  y3 <- morlet_wavelet(x, freqs, srate, wave_num, trend = trend, precision = "double")
  gc()
})


profile({
  y1 <- rave::wavelet(x, freqs, srate, wave_num, demean = TRUE)
})

# a <- y2[]
a <- y3$real[] + y3$imag[] * 1i
b <- t(y1$coef * exp(1i*y1$phase))

dif <- Mod(a-b)# / Mod(a)
max(dif)
idx <- which.max(dif)
a[idx]
b[idx]

# x <- rnorm(10000)
#
# list2env(list(freqs = 2:200, srate = 2000, wave_num = c(2,20), demean = TRUE), .GlobalEnv)
# ratio = (log(max(wave_num)) - log(min(wave_num))) / (log(max(freqs)) - log(min(freqs)))
# wavelet_cycles = round(exp((log(freqs) - log(min(freqs))) * ratio + log(min(wave_num))))
# d_l <- length(x)
#
# k1 <- ravetools:::wavelet_kernels(freqs, srate, wavelet_cycles)
# k2 <- rave:::wavelet_kernels(freqs, srate, wavelet_cycles)
# dev.off()
# sapply(1:length(k1$kernels), function(i){
#   max(Mod(k1$kernels[[i]] - k2[[i]]))
# })
# sapply(k1$kernels, length) - sapply(k2, length)
# length(k1$kernels[[2]])
# length(k2[[2]])
#
# wave_num <- wavelet_cycles
# k1 <- ravetools::wavelet_kernels2(freqs, srate, wave_num, d_l)
# rk <- function(data){
#   srate = round(srate);
#   if(length(wave_num) != length(freqs)){
#     # calculate wavelet cycles for each frequencies
#     ratio = (log(max(wave_num)) - log(min(wave_num))) / (log(max(freqs)) - log(min(freqs)))
#     wavelet_cycles = exp((log(freqs) - log(min(freqs))) * ratio + log(min(wave_num)))
#   }else{
#     wavelet_cycles = wave_num
#   }
#
#   # Instead of using fixed wave_cycles, use flex cycles
#   # lower num_cycle is good for low freq, higher num_cycle is good for high freq.
#   # wavelet_cycles = wave_num;
#   # lowest_freq = freqs[1];
#
#   f_l = length(freqs)
#   d_l = length(data)
#
#   # normalize data, and fft
#   if(demean){
#     fft_data = fftwtools::fftw_r2c(data - mean(data))
#   }else{
#     fft_data = fftwtools::fftw_r2c(data)
#   }
#
#
#   # wavelet window calc - each columns of final wave is a wavelet kernel (after fft)
#   # sts = wavelet_cycles / (2 * pi * freqs)
#   # wavelet_wins = cbind(-3 * sts, 3 * sts)
#
#   sapply(1:f_l, function(ii){
#     fq = freqs[ii]
#     cycles = wavelet_cycles[ii]
#     # standard error
#     st = cycles / (2 * pi * fq)
#
#     # calculate window size
#     wavelet_win = seq(-3 * st, 3 * st, by = 1/srate)
#
#     # half of window length
#     w_l_half = (length(wavelet_win) - 1) / 2
#
#     # wavelet 1: calc sinus in complex domain
#     tmp_sine = exp((0+1i) * 2 * pi * fq / srate * (-w_l_half:w_l_half))
#
#     # Gaussian normalization part
#     A = 1/sqrt(st*sqrt(pi))
#
#     # wavelet 2: calc gaussian wrappers
#     tmp_gaus_win = A * exp(-wavelet_win^2/(2 * (cycles/(2 * pi * fq))^2))
#
#     # wave kernel
#     tmp_wavelet = tmp_sine * tmp_gaus_win
#     tmp_wavelet
#
#     # padding
#     w_l = length(tmp_wavelet)
#     n_pre  = ceiling(d_l / 2) - floor(w_l/2)
#     n_post = d_l - n_pre - w_l
#     wvi = c(rep(0, n_pre), tmp_wavelet, rep(0, n_post))
#     fft_wave = Conj(fftwtools::fftw_c2c(wvi))
#
#     fft_wave
#   }) ->
#     fft_waves
#   fft_waves
# }
# k2 <- rk(x)
#
# range(Re(k1[] - k2))
# range(Im(k1[] - k2))
#
# y2 <- morlet_wavelet(x, freqs, srate, wave_num, demean)
# y1 <- rave::wavelet(x, freqs, srate, wave_num, demean)
#
# a <- y2[]
# b <- t(y1$coef * exp(1i*y1$phase))
#
# plot(a[2,], type = 'l'); lines(b[2,], col = 'red')
