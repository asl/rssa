AC_INIT([Rssa], 0.9.0)

dnl Select an optional include path, from a configure option
dnl or from an environment variable.
AC_ARG_WITH([fftw-include],
            AC_HELP_STRING([--with-fftw-include=INCLUDE_PATH],
                           [the location of FFTWv3 header files]),
            [fftw_include_path=$withval])

RSSA_CPPFLAGS="-I."
if test [ -n "$fftw_include_path" ] ; then
   RSSA_CPPFLAGS="-I. -I${fftw_include_path}"
else
   if test [ -n "${FFTW_INCLUDE}" ] ; then
      RSSA_CPPFLAGS="-I. -I${FFTW_INCLUDE}"
   fi
fi
     
dnl Ditto for a library path
AC_ARG_WITH([fftw-lib],
            AC_HELP_STRING([--with-fftw-lib=LIB_PATH],
                           [the location of FFTWv3 libraries]),
            [fftw_lib_path=$withval])
     
if test [ -n "$fftw_lib_path" ] ; then
   LIBS="-L$fftw_lib_path ${LIBS}"
else
   if test [ -n "${FFTW_LIBS}" ] ; then
      LIBS="-L${FFTW_LIBS} ${LIBS}"
   fi
fi

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

dnl Setup the compilers
AC_PROG_CC
AC_PROG_CPP

CPPFLAGS="${CPPFLAGS} ${RSSA_CPPFLAGS}"

dnl Check the headers can be found
AC_CHECK_HEADERS(fftw3.h,
                 AC_SUBST(HAVE_FFTW, 1),
                 AC_SUBST(HAVE_FFTW, 0))
     
AC_CHECK_LIB(fftw3, fftw_execute)

dnl Substitute RSSA_CPPFLAGS and LIBS
AC_SUBST(RSSA_CPPFLAGS)
AC_SUBST(LIBS)

dnl Do substitution in src/Makevars.in, src/config.h and R/init.R
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars R/init.R])
AC_OUTPUT
