
# Full Installation-from-Source Guide for `ravetools` 

The package `ravetools` contains `C++` code that requires general compilers such as `gcc` or `clang` to compile. In addition, the signal processing code requires [`FFTW3` library](https://www.fftw.org/), which could be installed easily. 

This guide contains three parts. In the first part, we will install proper compilers on your machine. If you have already installed them, please skip and proceed to the second part, in which you will install the `FFTW3` library. The last part simply installs `ravetools`.

#### Prerequisite

The development of `ravetools` is and will always be on the latest version of R. While I will try my best to maximize backward compatibility, it is unlikely for me to test `ravetools` on every single old versions of R. Choosing the latest R is always the best choice.


## Installing on Windows 

#### 1. Install Building Tools

* Please go to [this website](https://cran.r-project.org/bin/windows/Rtools/) and follow their instructions to download `Rtools`. As of `ravetools` is developed (R version `4.1.2`), the latest `Rtools` version is `4`. 

* Once `Rtools` is downloaded and installed, please open your R console, copy-paste the following R command and hit return:

```r
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)
```

#### 2. Install `FFTW3`


![In start menu, the `Rtools Bash` application icon is purple with a letter `M`](https://user-images.githubusercontent.com/216319/73364595-1fe28080-42ab-11ea-9858-ac8c660757d6.png)


* Open your Windows start menu, search `Rtools Bash`. If `Rtools` has been installed successfully, you will see an application with a purple icon and letter `M`. Open the application, you will see a terminal window.

* In the terminal, paste the following script and execute.

```
pacman -S  mingw-w64-{i686,x86_64}-fftw
```

* You will be prompted with a question `Proceed with installation? (Y/n)`. Please enter `Y` to agree. The `FFTW3` library will be added to your system.

#### 3. Install `ravetools`

* If you have previously opened any, please close all of them just in case. 
* Open your `R` or `RStudio`. In your newly opened `R` console, paste the following code and return line-by-line.

```r
if(system.file(package = 'remotes') == ""){ install.packages('remotes') }

remotes::install_github('dipterix/ravetools')
```


## Installing on MacOSX

#### 1. Check your R architecture

Open your R. In the start-up message, you will see text similar to these

```
R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: aarch64-apple-darwin20 (64-bit)
```

Alternatively, you can run the following R command:

```
R.version$platform
##> [1] "aarch64-apple-darwin20"
```

Check the architecture keywords after `Platform`. If you see `aarch64`, then you have ARM CPU (such as M1 chip), and an ARM version of R installed. If you see `x86_64`, then you have an Intel version of R installed (even if you have ARM chips, your R is still Intel-based).

> If you have ARM CPU (like M1-chip), it is always (strongly) recommended that you install R with the same architecture, because your compilers need to be consistent in architectures. If you have Intel-based R installed, then you have to set up a whole compiler chain that is Intel-based. (This tutorial does not cover such advanced case)

#### 2. Install command-line tools

* Open your `Terminal.app`. It can be easily found in `/Application` folder, or by pressing `command` key and the `space` key at the same time, type "terminal".

* Run `xcode-select --install` in the terminal and hit `return/enter` key. You will be prompted with an "Agreement" window. Please accept and wait the installation process to finish.

> If you see the message: `code-select: error: command line tools are already installed, use "Software Update" to install updates`, this means the command-line tools have been already installed on your machine in the past, and you are safe to proceed to the next step.

#### 3. Install `HomeBrew`

Please go to [HomeBrew's installation site](https://docs.brew.sh/Installation) to install `HomeBrew`. If you don't like to read, simply paste the following command into your terminal window and hit the `return` key:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

Please answer all the prompted questions and wait till `brew` to brew itself.

#### 4. Install `pkg-config` and `FFTW3`

Open a new terminal window. If you are using `ARM` CPU, then type

```
eval "$(/opt/homebrew/bin/brew shellenv)"
brew install pkg-config fftw
```

If you have Intel-based machine, type

```
eval "$(/usr/local/bin/brew shellenv)"
brew install pkg-config fftw
```

#### 5. Install `ravetools`

* If you have previously opened any, please close all of them just in case. 
* Open your `R` or `RStudio`. In your newly opened `R` console, paste the following code and return line-by-line.

```r
if(system.file(package = 'remotes') == ""){ install.packages('remotes') }

remotes::install_github('dipterix/ravetools')
```


## Installing on Ubuntu

* The following command will install essential build tools, system package configuration tools, and `FFTW3`

```
sudo apt-get install build-essential pkg-config libfftw3-dev
```

Please search their alternatives if you are using other Linux systems such as RedHat or CentOS.

* After installation, type `R` in the terminal, hit return/enter key to open R environment, and run the following R commands.

```r
if(system.file(package = 'remotes') == ""){ install.packages('remotes') }

remotes::install_github('dipterix/ravetools')
```

## Installing with Docker

Please check this [Dockerfile](https://github.com/dipterix/ravetools/blob/master/Dockerfile). It takes 24 minutes to compile everything from source with 2 CPUs, 1 GB RAM

