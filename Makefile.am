bin_PROGRAMS = 

lib_LTLIBRARIES=

libRVMnest_la_SOURCES = get_globalRVM_chi2.C get_RVM_chi2.C
libRVMnest_la_LIBADD = @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@

AM_CXXFLAGS = -std=gnu++11

if HAVE_MULTINEST
lib_LTLIBRARIES += libRVMnest.la
bin_PROGRAMS += globalRVMnest RVMnest precessnest precessnestmodel

libRVMnest_la_SOURCES += simpleParameters.C Parameters.C RVMnest_utils.C globalRVM_likelihood.C RVM_likelihood.C read_results.C precess_likelihood.C get_precessRVM_chi2.C  get_precessmRVM_chi2.C
libRVMnest_la_LIBADD += @FFTW_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

if HAVE_POLYCHORD
libRVMnest_la_SOURCES += globalRVM_likelihood_PC.C precess_likelihood_PC.C precessmodel_likelihood_PC.C polychord_c_interface.C
libRVMnest_la_LIBADD += @POLYCHORD_LIBS@
endif

if HAVE_MINUIT
lib_LTLIBRARIES += libRVMnuit.la
bin_PROGRAMS += globalRVMnuit RVMnuit

libRVMnuit_la_SOURCES = write_results.C naghi.c nrutil.c
libRVMnuit_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
libRVMnest_la_LIBADD += @MINUIT_LIBS@
endif

AM_CPPFLAGS = -I@top_srcdir@/ @FFTW_CFLAGS@  @CFITSIO_CFLAGS@ @PSRCHIVE_CFLAGS@  @GLIB_CFLAGS@ @MINUIT_CFLAGS@ @MULTINEST_CFLAGS@ 
if HAVE_POLYCHORD
AM_CPPFLAGS += @POLYCHORD_CFLAGS@
endif


globalRVMnest_SOURCES = globalRVMnest.C
globalRVMnest_LDADD = libRVMnest.la @PSRCHIVE_LIBS@

RVMnest_SOURCES = RVMnest.C
RVMnest_LDADD = libRVMnest.la

precessnest_SOURCES = precessnest.C
precessnest_LDADD = libRVMnest.la @PSRCHIVE_LIBS@

precessnestmodel_SOURCES = precessnestmodel.C
precessnestmodel_LDADD = libRVMnest.la @PSRCHIVE_LIBS@

if HAVE_POLYCHORD
globalRVMnest_LDADD +=  @POLYCHORD_LIBS@
endif

if HAVE_MINUIT
globalRVMnuit_SOURCES = globalRVMnuit.C
globalRVMnuit_LDADD = libRVMnuit.la

RVMnuit_SOURCES = RVMnuit.C
RVMnuit_LDADD = libRVMnuit.la
endif
