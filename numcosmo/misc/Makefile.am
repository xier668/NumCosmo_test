## Process this file with automake to produce Makefile.in

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS) \
	 $(SUNDIALS_CFLAGS) \
	 $(GSL_CFLAGS) \
	 $(FFTW3_CFLAGS) \
	 $(CFITSIO_CFLAGS) \
	 $(NLOPT_CFLAGS) \
	 -I$(top_srcdir)/numcosmo \
	 -I$(builddir)/.. \
	 $(VISIBILITY_HIDDEN_FLAG)

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

AM_LDFLAGS = $(VISIBILITY_HIDDEN_FLAG)

SUBDIRS = .

noinst_LTLIBRARIES = libnumcosmo_misc.la

libnumcosmo_misc_la_LIBADD = \
	$(SUNDIALS_LIBS) \
	$(FFTW3_LIBS) \
	$(GLIB_LIBS) \
	$(NLOPT_LIBS) \
	$(GSL_LIBS) \
	$(CFITSIO_LIBS)

libnumcosmo_misc_la_SOURCES = \
	Faddeeva.c      \
	Faddeeva.h      \
	rcm.c           \
	rcm.h           \
	LowRankQP.c     \
	libqp_gsmo.c    \
	libqp_splx.c    \
	libqp.h         \
	nnls.c          \
	cqpminimizer.c  \
	initial_point.c \
	mg_pdip.c       \
	gsl_cqp.h       \
 	CNearTree.c	\
	CNearTree.h	\
	CVector.h	\
	CVector.c       \
	rhrand.h        \
	initial_point.h

EXTRA_DIST =

CLEANFILES = 
