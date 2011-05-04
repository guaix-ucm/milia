/*
 * Copyright 2008-2011 Sergio Pascual
 *
 * This file is part of Milia
 *
 * Milia is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milia is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milia.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// $Id$

#ifndef MILIA_FLRW_TEST_H
#define MILIA_FLRW_TEST_H

#include "milia/exception.h"
#include "FlrwTestDataMixin.h"

#include <cppunit/extensions/HelperMacros.h>

class FlrwTest : public CppUnit::TestFixture, public FlrwTestDataMixin
{
    CPPUNIT_TEST_SUITE(FlrwTest);
    CPPUNIT_TEST_EXCEPTION(testHubbleZeroThrows, milia::exception);
    CPPUNIT_TEST_EXCEPTION(testHubbleLessThanZeroThrows, milia::exception);
    CPPUNIT_TEST(testLuminosityDistance);
    CPPUNIT_TEST(testAngularDistance);
    CPPUNIT_TEST(testComovingTransverseDistance);
    CPPUNIT_TEST(testComovingDistance);
    CPPUNIT_TEST(testAge);
    CPPUNIT_TEST(testComovingVolume);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();

    void tearDown();

    /** Tests hubble = 0 */
    void testHubbleZeroThrows();

    /** Tests hubble < 0 */
    void testHubbleLessThanZeroThrows();

    void testLuminosityDistance();

    void testComovingDistance();

    void testComovingTransverseDistance();

    void testAngularDistance();

    void testAge();

    void testComovingVolume();

private:
    static const double ms_rel_tol = 1e-4;
    static const double ms_z = 0.1;
    static const double ms_hubble = 65;
};


#endif // MILIA_FLRW_TEST_H
