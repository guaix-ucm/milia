## $Id$

lib_LTLIBRARIES = libmilia.la
libmilia_la_SOURCES = flrw.cc flrw_cached.cc flrw_distance.cc flrw_age.cc exception.cc \
    schechter.cc

## Version info is in current : revision : age form
## A library supports interfaces from current downto current - age
## Revision is the version of the current interface

## If the source code has changed, rev += 1

libmilia_la_LDFLAGS = -version-info 1:0:1
libmilia_la_LIBADD = $(GSL_LIBS)

pkginclude_HEADERS = metric.h flrw.h flrw_cached.h exception.h schechter.h

AM_CPPFLAGS = $(GSL_CFLAGS) -I$(top_srcdir)

