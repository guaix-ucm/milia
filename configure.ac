# Process this file with autoconf to produce a configure script.

# $Id$

AC_INIT(milia,0.1.0)
AM_INIT_AUTOMAKE([foreign])

AC_REVISION($Revision$)

AC_PROG_CXX
AC_PROG_LIBTOOL

PKG_CHECK_MODULES([STUFF], [gsl])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_FILES([Makefile 
		milia.pc
		milia/Makefile
		tests/Makefile])
AC_OUTPUT
