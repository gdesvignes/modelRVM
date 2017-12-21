bin_PROGRAMS = 

lib_LTLIBRARIES=

libRVMnest_la_SOURCES = get_globalRVM_chi2.C
libRVMnest_la_LIBADD = @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@

if HAVE_MULTINEST
lib_LTLIBRARIES += libRVMnest.la
bin_PROGRAMS += globalRVMnest RVMnest

libRVMnest_la_SOURCES += simpleParameters.C Parameters.C RVMnest_utils.C globalRVM_likelihood.C RVM_likelihood.C read_results.C
libRVMnest_la_LIBADD += @FFTW_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

if HAVE_POLYCHORD
libRVMnest_la_SOURCES += globalRVM_likelihood_PC.C polychord_c_interface.C
libRVMnest_la_LIBADD += @POLYCHORD_LIBS@
endif

if HAVE_MINUIT
lib_LTLIBRARIES += libRVMnuit.la
bin_PROGRAMS += globalRVMnuit RVMnuit

libRVMnuit_la_SOURCES = write_results.C naghi.c nrutil.c
libRVMnuit_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

AM_CPPFLAGS = -I@top_srcdir@/ @FFTW_CFLAGS@  @CFITSIO_CFLAGS@ @PSRCHIVE_CFLAGS@  @GLIB_CFLAGS@ @MINUIT_CFLAGS@ @MULTINEST_CFLAGS@ 
if HAVE_POLYCHORD
AM_CPPFLAGS += @POLYCHORD_CFLAGS@
endif


if HAVE_MULTINEST
globalRVMnest_SOURCES = globalRVMnest.C
globalRVMnest_LDADD = libRVMnest.la

RVMnest_SOURCES = RVMnest.C
RVMnest_LDADD = libRVMnest.la
endif

if HAVE_MINUIT
globalRVMnuit_SOURCES = globalRVMnuit.C
globalRVMnuit_LDADD = libRVMnuit.la

RVMnuit_SOURCES = RVMnuit.C
RVMnuit_LDADD = libRVMnuit.la
endif
