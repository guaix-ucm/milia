if TESTS_ENABLED
TESTS = flrw_test

check_PROGRAMS = flrw_test

flrw_test_SOURCES = flrw_test.cc FlrwTest.h FlrwTest.cc FlrwTestData.cc \
  FlrwTestData.h FlrwAge.h FlrwAge.cc

AM_CPPFLAGS = -I$(top_srcdir) -I$(CPPUNIT_CFLAGS)
LDADD = $(top_builddir)/milia/libmilia.la $(CPPUNIT_LIBS)
else
TESTS = no_tests.sh
XFAIL_TESTS = no_tests.sh

dist_check_SCRIPTS = no_tests.sh
endif
