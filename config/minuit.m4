#
# SWIN_LIB_MINUIT([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the package MultiNest
# by F. Feroz, M.P. hobson, and M. Bridges
#
# MINUIT_CFLAGS - autoconfig variable with flags required for compiling
# MINUIT_LIBS   - autoconfig variable with flags required for linking
# HAVE_MINUIT   - automake conditional
# HAVE_MINUIT   - pre-processor macro in config.h
#
# This macro tries to link a test program, by using something like
#
#    -L$MINUIT_HOME -lnest3 -llapack -blas
#
# Notice that the environment variable MINUIT_HOME is used. In the case the
# library 'libnest' is not installed in a default location, let MINUIT_HOME
# point to the location of libnest
#
#  MINUIT_LIBS="$MINUIT_LIBS -lnest3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"
#  MINUIT_LIBS="$MINUIT_LIBS -lnest3 -llapack"
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_MINUIT],
[
  AC_PROVIDE([SWIN_LIB_MINUIT])

  SWIN_PACKAGE_OPTIONS([minuit])

  if test x"$MINUIT_HOME" != x; then
    MINUIT_CFLAGS="-I$MINUIT_HOME/include"
    MINUIT_LIBS="-L$MINUIT_HOME/lib"
    have_minuit=yes
  else 
    MINUIT_CFLAGS="-I$with_minuit_include_dir"
    MINUIT_LIBS="-L$with_minuit_lib_dir"
  fi
  ac_save_CFLAGS="$CFLAGS"
  ac_save_CXXFLAGS="$CXXFLAGS"
  ac_save_LIBS="$LIBS"

  MINUIT_LIBS="$MINUIT_LIBS -lMinuit2 -fopenmp"

  LIBS="$ac_save_LIBS $MINUIT_LIBS"
  CFLAGS="$ac_save_CFLAGS $MINUIT_CFLAGS"
  CXXFLAGS="$ac_save_CXXFLAGS $MINUIT_CFLAGS"


  AC_MSG_CHECKING([for Minuit installation])

  AC_LANG_PUSH(C++)

  

      AC_TRY_RUN(
          [


	#include "Minuit2/FunctionMinimum.h"
	#include "Minuit2/MnPrint.h"
	#include "Minuit2/VariableMetricMinimizer.h"
	#include "Minuit2/MnMigrad.h"
	#include "Minuit2/MnMinos.h"

	using namespace ROOT::Minuit2;

	int main() {
	      MnUserParameters upar;
	      upar.Add("x", 1., 0.1);
	  return 0;
	    } /* main */
		  ],
		  [
		      have_minuit=yes
		  ],
		  [
		  ])
      AC_LANG_POP(C++)

  AC_MSG_RESULT($have_minuit)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"
  CXXFLAGS="$ac_save_CXXFLAGS"

  if test x"$have_minuit" = xyes; then
    AC_DEFINE([HAVE_MINUIT], [1], [Define to 1 if you have the MINUIT library])
    [$1]
  else
    AC_MSG_WARN([Minuit code will not be compiled. This only affects MultiNest plugins.])

    if test x"$have_minuit" = xyes; then
      AC_MSG_WARN([***************************************************************])
      AC_MSG_WARN([Minuit code can be compiled, but it cannot be run. ])
      AC_MSG_WARN([Please set (DY)LD_LIBRARY_PATH correctly.])
      AC_MSG_WARN([***************************************************************])
    else
      if test x"$MINUIT_HOME" = x; then
        AC_MSG_NOTICE([Minuit plugins will not be compiled. If you need this, set MINUIT_HOME.])
# Reduced warning message size
#        AC_MSG_WARN([***************************************************************])
#        AC_MSG_WARN([Please set the MINUIT_HOME environment variable, or set the LIBRARY_PATH variable])
#        AC_MSG_WARN([       Note: LIBRARY_PATH is used at linking stage, (DY)LD_LIBRARY_PATH is used at runtime])
#        AC_MSG_WARN([***************************************************************])
      fi
    fi

    MINUIT_CFLAGS=""
    MINUIT_LIBS=""
    [$2]
  fi

  AC_SUBST(MINUIT_CFLAGS)
  AC_SUBST(MINUIT_LIBS)
  AM_CONDITIONAL(HAVE_MINUIT, [test x"$have_minuit" = xyes])
])

