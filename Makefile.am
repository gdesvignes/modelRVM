bin_PROGRAMS = 

lib_LTLIBRARIES=libRVMnest.la libRVMnuit.la


if HAVE_MULTINEST
bin_PROGRAMS += globalRVMnest RVMnest

libRVMnest_la_SOURCES = Parameters.C RVMnest_utils.C globalRVM_likelihood.C RVM_likelihood.C
libRVMnest_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

if HAVE_MINUIT
bin_PROGRAMS += RVMnuit

libRVMnuit_la_SOURCES = write_results.C naghi.c nrutil.c
libRVMnuit_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

AM_CPPFLAGS = -I@top_srcdir@/ @FFTW_CFLAGS@  @CFITSIO_CFLAGS@ @PSRCHIVE_CFLAGS@  @GLIB_CFLAGS@ @MINUIT_CFLAGS@ @MULTINEST_CFLAGS@ 

if HAVE_MULTINEST
globalRVMnest_SOURCES = globalRVMnest.C
globalRVMnest_LDADD = libRVMnest.la

RVMnest_SOURCES = RVMnest.C
RVMnest_LDADD = libRVMnest.la
endif

if HAVE_MINUIT
RVMnuit_SOURCES = RVMnuit.C
RVMnuit_LDADD = libRVMnuit.la
endif
