# install.packages(c("fftw", "microbenchmark", "fftwtools", "lobstr"))

# The performance has a phase-transition when
# the signal length is around 100000
x <- rnorm(100000) + 1; lobstr::obj_size(x); x[[1]]


opt <- 0L
p <- fftw::planFFT(length(x), opt)
microbenchmark::microbenchmark({
  fftw::FFT(x, plan = p)
},{
  fftw::FFT(x)
}, {
  stats::fft(x)
}, {
  ravetools:::fftw_r2c(x, inplace = FALSE, fftwplanopt = opt)
}, {
  ravetools:::fftw_r2c(x, inplace = TRUE, fftwplanopt = opt)
}, {
  fftwtools::fftw_r2c(x)
}, times = 10, check = function(v){
  y <- v[[1]]
  all(sapply(v, function(x){
    dif <- max(Mod(x-y))
    dif < 1e-7
  }))
}); x[[1]]

