# Process this file with autoconf to produce a configure script.

# $Id$

AC_INIT(milia,0.1.0)
AM_INIT_AUTOMAKE([foreign])

AC_REVISION($Revision$)

AC_PROG_CXX
AC_PROG_LIBTOOL

PKG_CHECK_MODULES([STUFF], [gsl])

# For the tests
PKG_CHECK_MODULES([CPPUNIT], [cppunit], [testen=1], [
AC_MSG_WARN([Cppunit is needed by the tests. Checking is disabled])
testen=0])
AM_CONDITIONAL([TESTS_ENABLED], [test x$testen = x1])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_FILES([Makefile 
		milia.pc
		milia/Makefile
		examples/Makefile
		tests/Makefile])
AC_OUTPUT
