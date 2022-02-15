
AC_PREREQ([2.61])

## Update this version number along with DESCRIPTION
AC_INIT
AC_CONFIG_SRCDIR([ravetools])

# Find R home
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
AC_REQUIRE_CPP
AC_PROG_CC

# RCPP_PKG_LIBS=`${R_HOME}/bin/Rscript -e "Rcpp:::LdFlags()"`

if test `uname` = "Darwin" ; then
  darwin="yes"
else
  darwin="no"
fi

## Check pkg-config
AC_DEFUN([AC_PROG_PKGCONFIG], [AC_CHECK_PROG(PKGCONFIG,pkg-config,yes)])
AC_PROG_PKGCONFIG

if test x"${PKGCONFIG}" != x"yes"; then
   AC_MSG_ERROR([pkg-config not found.])
fi

## Flags
fftw3_cflags=`pkg-config --cflags fftw3`
fftw3_libs=`pkg-config --libs fftw3`


## --- Output ----------------------------------------------------------------
## Substitute and generate Makevars
AC_SUBST([PKG_CFLAGS],["${fftw3_cflags} ${PKG_CFLAGS} -I../inst/include"])
#  ${RCPP_PKG_LIBS}
AC_SUBST([PKG_LIBS],["${fftw3_libs} ${PKG_LIBS}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "
================================
compiler flags: ${PKG_CFLAGS}
link argument:  ${PKG_LIBS}
================================
"
