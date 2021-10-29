
time <- seq(0, 3, by = 0.001)
x <- sin(time * 20*pi) + exp(-time^2) * cos(time * 10*pi)

res <- multitaper(
  data = x, fs = 1000, frequency_range = c(0,15),
  time_bandwidth=1,
  window_params=c(0.5,0.01), num_tapers = 1, nfft = 500
)
res$frequency
wave <- morlet_wavelet(data = x, freqs = res$frequency, srate = 1000, wave_num = c(5,8), demean = TRUE)

image(
  x = res$time,
  y = res$frequency,
  z = res$spec,
  # z = 10 * log10(res$spec),
  xlab = "Time (s)",
  ylab = 'Frequency (Hz)',
  col = matlab_palette()
)
image(Mod(wave[seq(1, 3001, by = 10),])^2)
plot(rowMeans(Mod(wave[seq(1, 3001, by = 10),])^2), type='l')

fields::image.plot(
  x = res$time,
  y = res$frequency,
  # res$spec,
  10 * log10(res$spec),
  xlab = "Time (s)",
  ylab = 'Frequency (Hz)'
)


multitaper_spectrogram_R(
  x, 100, frequency_range = c(0,15),
  time_bandwidth=1.5,
  window_params=c(1.3,0.01)
)




data <- sin(100 * seq(0, 10, length.out = 10000)) + rnorm(10000)

win_size <- 0.1
max_freq <- 40
fs <- 1000
freq_steps <- 2
time_bandwidth <- win_size * max_freq / 2
conf <- multitaper_config(
  length(data), 1000, frequency_range=c(0,30),
  time_bandwidth=2,
  window_params=c(1,0.01)
)
conf$nfft

# nfft <- fs / freq_steps
res <- multitaper(
  data, 1000, frequency_range=c(0,30),
  time_bandwidth=2,
  nfft = NA,
  window_params=c(5,0.01)
)

fields::image.plot(
  x = res$time,
  y = res$frequency,
  # res$spec,
  10 * log10(res$spec),
  xlab = "Time (s)",
  ylab = 'Frequency (Hz)'
)


multitaper_spectrogram_R(
  data,
  1000,
  frequency_range=c(0,max_freq),
  time_bandwidth = 1,
  window_params=c(1,0.02)
)

