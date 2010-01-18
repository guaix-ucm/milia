## $Id$

lib_LTLIBRARIES = libmilia.la
libmilia_la_SOURCES = flrw_nat.cc flrw.cc flrw_nat_distance.cc flrw_nat_age.cc \
    exception.cc flrw_prec.h schechter.cc

## Version info is in current : revision : age form
## A library supports interfaces from current downto current - age
## Revision is the version of the current interface

## If the source code has changed, rev += 1

libmilia_la_LDFLAGS = -version-info 3:0:0
libmilia_la_LIBADD = $(STUFF_LIBS)

pkginclude_HEADERS = metric.h flrw.h flrw_nat.h exception.h schechter.h

AM_CPPFLAGS = -I$(top_srcdir) $(STUFF_CFLAGS) 

