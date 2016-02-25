bin_PROGRAMS = globalRVMnest

lib_LTLIBRARIES=libRVMnest.la


c_sources = Parameters.C RVMnest_utils.C globalRVM_likelihood.C

libRVMnest_la_SOURCES = $(c_sources)

libRVMnest_la_LIBADD = @FFTW_LIBS@ @CFITSIO_LIBS@ @PSRCHIVE_LIBS@ @GLIB_LIBS@ @MINUIT_LIBS@

AM_CPPFLAGS = -I@top_srcdir@/ @FFTW_CFLAGS@  @CFITSIO_CFLAGS@ @PSRCHIVE_CFLAGS@  @GLIB_CFLAGS@ @MINUIT_CFLAGS@

globalRVMnest_SOURCES = globalRVMnest.C
globalRVMnest_LDADD = libRVMnest.la
