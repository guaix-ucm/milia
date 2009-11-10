# Process this file with autoconf to produce a configure script.

# $Id$

AC_INIT([milia],[0.3.0])
AM_INIT_AUTOMAKE([foreign])
LT_INIT([disable-shared])

AC_REVISION($Revision$)

AC_PROG_CXX

AC_LANG([C++])
AC_CHECK_HEADER([boost/function.hpp],[],[AC_MSG_ERROR([no boost])])
PKG_CHECK_MODULES([GSL], [gsl])
# For the tests
PKG_CHECK_MODULES([CPPUNIT], [cppunit], [testen=1], [
AC_MSG_WARN([Cppunit is needed by the tests. Checking is disabled])
testen=0])
AM_CONDITIONAL([TESTS_ENABLED], [test x$testen = x1])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_FILES([Makefile 
         	 milia.pc
         	 Doxyfile
		 milia/Makefile
		 examples/Makefile
		 tests/Makefile])
AC_OUTPUT
