dat <- raveio::read_mat('~/rave_data/raw_dir/YAB/008/YABDatafile008_ch14.mat')
dat <- raveio::read_mat('/Volumes/BeauchampServe/rave_data/raw/YAB/008/YABDatafile008_ch14.mat')
data <- as.vector(dat$analogTraces)[10000 + 0:10000]
# data <- rnorm(1000)
win_size <- 0.05
max_freq <- 200
fs <- 2000
freq_steps <- 2
time_bandwidth <- win_size * max_freq / 2
conf <- multitaper_config(
  length(data), fs, frequency_range=c(0,max_freq),
  time_bandwidth=time_bandwidth,
  num_tapers=NULL, window_params=c(win_size,0.01),
  nfft=fs / freq_steps, detrend_opt='linear'
)
conf

res <- multitaper(data, fs, frequency_range=c(0,max_freq),
           time_bandwidth=time_bandwidth,
           num_tapers=NULL, window_params=c(win_size,0.01),
           nfft=fs / freq_steps, detrend_opt='linear', plot_on=TRUE, verbose=TRUE)


fields::image.plot(
  x = res$time,
  y = res$frequency,
  # res$spec,
  10 * log10(res$spec),
  xlab = "Time (s)",
  ylab = 'Frequency (Hz)'
)

plot(colMeans(log(res$spec)), type = 'l')


