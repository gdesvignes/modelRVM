bin_PROGRAMS = 

lib_LTLIBRARIES=libRVMnest.la libRVM.la


if HAVE_MULTINEST
bin_PROGRAMS += globalRVMnest

libRVMnest_la_SOURCES = Parameters.C RVMnest_utils.C globalRVM_likelihood.C
libRVMnest_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

if HAVE_MINUIT
bin_PROGRAMS += fitRVMKJ

libRVM_la_SOURCES = write_results.C naghi.c nrutil.c
libRVM_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@ @MULTINEST_LIBS@
endif

AM_CPPFLAGS = -I@top_srcdir@/ @FFTW_CFLAGS@  @CFITSIO_CFLAGS@ @PSRCHIVE_CFLAGS@  @GLIB_CFLAGS@ @MINUIT_CFLAGS@ @MULTINEST_CFLAGS@ 

if HAVE_MULTINEST
globalRVMnest_SOURCES = globalRVMnest.C
globalRVMnest_LDADD = libRVMnest.la
endif

if HAVE_MINUIT
fitRVMKJ_SOURCES = fit_RVMKJ.C
fitRVMKJ_LDADD = libRVM.la
endif
