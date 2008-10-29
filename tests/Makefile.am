## $Id$

TESTS = main_test

check_PROGRAMS = main_test

main_test_SOURCES = main_test.cc FlrwTest.h FlrwTest.cc FlrwTestData.cc

AM_CPPFLAGS = -I$(top_srcdir) -I$(CPPUNIT_CFLAGS)
LDADD = $(top_builddir)/milia/libmilia.la $(STUFF_LIBS) $(CPPUNIT_LIBS)
