lib_LTLIBRARIES = libmilia.la
libmilia_la_SOURCES = flrw.cc flrw_prec.h metric.cc\
    flrw_nat.cc flrw_nat_distance.cc flrw_nat_age.cc util.cc util.h

    

## Version info is in current : revision : age form
## A library supports interfaces from current downto current - age
## Revision is the version of the current interface

## If the source code has changed, rev += 1

libmilia_la_LDFLAGS = -version-info 4:0:0
libmilia_la_LIBADD = $(GSL_LIBS) $(BOOST_LDFLAGS)

pkginclude_HEADERS = metric.h flrw.h flrw_nat.h 

AM_CPPFLAGS = $(GSL_CFLAGS) $(BOOST_CPPFLAGS) -I$(top_srcdir)

