#
# SWIN_LIB_POLYCHORD([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
#
# This m4 macro checks availability of the package PolyChord
# by Will Handley, Mike Hobson & Anthony Lasenby
#
# POLYCHORD_CFLAGS - autoconfig variable with flags required for compiling
# POLYCHORD_LIBS   - autoconfig variable with flags required for linking
# HAVE_POLYCHORD   - automake conditional
# HAVE_POLYCHORD   - pre-processor macro in config.h
#
# This macro tries to link a test program, by using something like
#
#    -L$POLYCHORD_DIR -lchord
#
# Notice that the environment variable POLYCHORD_DIR is used. In the case the
# library 'libnest' is not installed in a default location, let POLYCHORD_DIR
# point to the location of libnest
#
#  POLYCHORD_LIBS="$POLYCHORD_LIBS -lnest3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread"
#  POLYCHORD_LIBS="$POLYCHORD_LIBS -lnest3 -llapack"
# ----------------------------------------------------------
AC_DEFUN([SWIN_LIB_POLYCHORD],
[
  AC_PROVIDE([SWIN_LIB_POLYCHORD])
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([ACX_LAPACK])
  AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])

  POLYCHORD_CFLAGS=""
  POLYCHORD_LIBS=""

  if test x"$POLYCHORD_DIR" != x; then
    POLYCHORD_CFLAGS="-I$POLYCHORD_DIR"
    POLYCHORD_LIBS="-L$POLYCHORD_DIR"
  fi

  POLYCHORD_LIBS="$POLYCHORD_LIBS $BLAS_LIBS $LIBS_LAPACK -lchord"

  ac_save_CFLAGS="$CFLAGS"
  ac_save_LIBS="$LIBS"
  LIBS="$ac_save_LIBS $POLYCHORD_LIBS"
  CFLAGS="$ac_save_CFLAGS $POLYCHORD_CFLAGS"

  AC_CHECK_LIB(chord, polychord_c_interface,
            [

                compile_polychord=yes
            ],
            [
                AC_CHECK_LIB(chord, polychord_c_interface,
                    [
                        compile_polychord=yes
                    ])
            ])

  if test x"$compile_polychord" = xyes; then
      AC_MSG_NOTICE(Running PolyChord to confirm runtime linking...)
      AC_LANG_PUSH(C++)
      AC_TRY_RUN(
          [

#include <stdio.h>
#include <string.h>
#include <string>

extern "C" void polychord_c_interface( double (*)(double*,int,double*,int,void*), void (*)(double*,double*,int,void*), int, int, bool, int, double, i\
nt, double, bool, bool, bool, bool, bool, bool, bool, bool, bool, int, int, int, char*, char*,int,double*,int*);

struct Settings
{
    int nDims;
    int nDerived;
    int nlive;
    int num_repeats;
    bool do_clustering;
    int feedback;
    double precision_criterion;
    int max_ndead;
    double boost_posterior;
    bool posteriors;
    bool equals;
    bool cluster_posteriors;
    bool write_resume;
    bool write_paramnames;
    bool read_resume;
    bool write_stats;
    bool write_live;
    bool write_dead;
    int update_files;
    std::string base_dir;
    std::string file_root;
    //Settings(int _nDims=0,int _nDerived=0);
};

void run_polychord( 
		   double (*c_loglikelihood_ptr)(double*,int,double*,int,void*), 
		   void (*c_prior_ptr)(double*,double*,int,void*), 
        Settings s)
{
    // Ridiculous gubbins for passing strings between C and FORTRAN
    char * base_dir = new char[s.base_dir.size()+1];
    std::copy(s.base_dir.begin(),s.base_dir.end(),base_dir);
    base_dir[s.base_dir.size()] = '\0';

    char * file_root = new char[s.file_root.size()+1];
    std::copy(s.file_root.begin(),s.file_root.end(),file_root);
    file_root[s.file_root.size()] = '\0';

    int ngrade = 1;
    double grade_frac[ngrade];
    grade_frac[0] = 1.;
    int grade_dims[ngrade];
    grade_dims[0] = s.nDims;

    polychord_c_interface( 
            c_loglikelihood_ptr, 
            c_prior_ptr, 
            s.nlive, 
            s.num_repeats,
            s.do_clustering,
            s.feedback,
            s.precision_criterion,
            s.max_ndead,
            s.boost_posterior,
            s.posteriors,
            s.equals,
            s.cluster_posteriors,
            s.write_resume,
            s.write_paramnames,
            s.read_resume,
            s.write_stats,
            s.write_live,
            s.write_dead,
            s.update_files,
            s.nDims,
            s.nDerived,
            base_dir,
            file_root,
            ngrade,
            grade_frac,
            grade_dims);

    delete[] base_dir;
    delete[] file_root;
}

	void run_polychord( double (*loglikelihood)(double*,int,double*,int,void *), void (*prior)(double*,double*,int,void *), Settings);

	double LogLike(double theta[], int nDims, double phi[], int nDerived, void *context) {
	  return -5 * (theta[0]*theta[0] + theta[1]*theta[1]);	 
	}

	void prior (double cube[], double theta[], int nDims,  void *context) {
	  theta[0] = cube[0];
	  theta[1] = cube[1];
	}


	int main(int argc, const char **argv) {
	  Settings settings;
          settings.nDims         = 2;
          settings.nDerived      = 1;
          settings.nlive       	 = 10;
          settings.num_repeats   = settings.nDims*5;
          settings.do_clustering = false;
          settings.precision_criterion = 0.1;
          settings.base_dir.assign("./");
          settings.file_root     = "test";
          settings.write_resume  = false;
          settings.read_resume   = false;
          settings.write_live    = false;
          settings.write_dead    = false;
          settings.write_stats   = false;
          settings.equals        = false;
          settings.posteriors    = false;
          settings.cluster_posteriors = false;
          settings.feedback      = 2;
          settings.update_files  = settings.nlive;
          settings.boost_posterior= 5.0;

	run_polychord(LogLike, prior, settings) ;


        return 0;
    } /* main */
          ],
          [
              AC_MSG_NOTICE(Running PolyChord successful)
              have_polychord=yes
          ],
          [
              AC_MSG_NOTICE(Running PolyChord failed)
          ])
      AC_LANG_POP(C++)
  fi

  AC_MSG_CHECKING([for PolyChord installation])
  AC_MSG_RESULT($have_polychord)

  LIBS="$ac_save_LIBS"
  CFLAGS="$ac_save_CFLAGS"

  if test x"$have_polychord" = xyes; then
    AC_DEFINE([HAVE_POLYCHORD], [1], [Define to 1 if you have the POLYCHORD library])
    [$1]
  else
    AC_MSG_WARN([PolyChord code will not be compiled. This only affects PolyChord plugins.])

    if test x"$compile_polychord" = xyes; then
      AC_MSG_WARN([***************************************************************])
      AC_MSG_WARN([PolyChord code can be compiled, but it cannot be run. This most])
      AC_MSG_WARN([likely means that the PolyChord library libnest3.so is not in the])
      AC_MSG_WARN([path. Please set (DY)LD_LIBRARY_PATH correctly.])
      AC_MSG_WARN([***************************************************************])
    else
      if test x"$POLYCHORD_DIR" = x; then
        AC_MSG_NOTICE([PolyChord plugins will not be compiled. If you need this, set POLYCHORD_DIR.])
# Reduced warning message size
#        AC_MSG_WARN([***************************************************************])
#        AC_MSG_WARN([Please set the POLYCHORD_DIR environment variable, or set the LIBRARY_PATH variable])
#        AC_MSG_WARN([       Note: LIBRARY_PATH is used at linking stage, (DY)LD_LIBRARY_PATH is used at runtime])
#        AC_MSG_WARN([***************************************************************])
      fi
    fi

    POLYCHORD_CFLAGS=""
    POLYCHORD_LIBS=""
    [$2]
  fi

  AC_SUBST(POLYCHORD_CFLAGS)
  AC_SUBST(POLYCHORD_LIBS)
  AM_CONDITIONAL(HAVE_POLYCHORD, [test x"$have_polychord" = xyes])
])

