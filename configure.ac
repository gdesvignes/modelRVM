AC_INIT([RVMnest], [1.0], [gdesvignes@mpifr-bonn.mpg.de])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall foreign])

AM_PROG_AR


AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77
AC_PROG_FC
AC_PROG_F77_C_O

AC_PROG_INSTALL
AC_PROG_LIBTOOL

AC_CHECK_LIB([m],[cos])


# Checks for essential libraries


ACX_BLAS(have_blas=yes,have_blas=no)
AC_MSG_CHECKING([for BLAS])
AC_MSG_RESULT($have_blas)
AM_CONDITIONAL(HAVE_BLAS, [test x$have_blas = xyes])

ACX_LAPACK(have_lapack=yes, have_lapack=no)
AC_MSG_CHECKING([for LAPACK])
AC_MSG_RESULT($have_lapack)
AM_CONDITIONAL(HAVE_LAPACK, [test x"$have_lapack" = xyes])

AX_MPI(have_mpi=yes,have_mpi=no)
AC_MSG_CHECKING([for MPI])
AC_MSG_RESULT($have_mpi)
AM_CONDITIONAL(HAVE_MPI, [test x"$have_mpi" = xyes])

PKG_CHECK_MODULES([GLIB], glib-2.0 >= 2.16)
#
AC_MSG_NOTICE([Adding glib-2.0 CFLAGS ${GLIB_CFLAGS}])
AM_CFLAGS="${GLIB_CFLAGS} ${AM_CFLAGS}"
AC_MSG_NOTICE([Adding glib-2.0 LIBS ${GLIB_LIBS}])
LIBS="${GLIB_LIBS} ${LIBS}"

SWIN_LIB_MINUIT
SWIN_LIB_FFTW
SWIN_LIB_MKL
SWIN_LIB_CFITSIO
SWIN_LIB_PSRCHIVE([],AC_MSG_ERROR([PSRCHIVE Library not found.
Please see http://psrchive.sourceforge.net
]))
SWIN_LIB_MULTINEST
SWIN_LIB_POLYCHORD



AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
])
AC_OUTPUT

