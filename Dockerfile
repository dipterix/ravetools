FROM rocker/r-base:latest

LABEL maintainer="Zhengjia Wang <dipterix.wang@gmail.com>" \
  org.label-schema.license="GPL-3.0" \
  org.label-schema.vcs-url="https://github.com/dipterix/ravetools" \
  org.label-schema.vendor="Dipterix"

# Install system libraries...
# If you only want to build `ravetools`, some of them are required
# However, if you plan to use RAVE, please don't miss any

RUN apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    build-essential \
    file \
    git \
    libsodium-dev \
    libffi-dev \
    libbz2-dev \
    libpcre2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libhdf5-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libpng-dev \
    libssl-dev \
    libssh2-1-dev \
    libtiff5-dev \
    libv8-dev \
    libxml2-dev \
    psmisc \
    procps \
    pkg-config \
    sudo \
    wget \
    zlib1g-dev \
  && install2.r --error --skipinstalled --deps TRUE \
    signal Rcpp RcppParallel remotes digest waveslim filearray docopt \
  # Get NCPUs
  && ncpus=$(Rscript --no-save -e "cat(parallel::detectCores())") \
  && echo "Installing Github version of ravetools with $ncpus CPUs" \
  && Rscript --verbose --no-save -e "remotes::install_github('dipterix/ravetools', upgrade = FALSE, force = TRUE, Ncpus = $ncpus)"

# Local build only
# docker build --tag ravetools:devel .


