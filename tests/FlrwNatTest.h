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

#ifndef MILIA_FLRW_NAT_TEST_H
#define MILIA_FLRW_NAT_TEST_H

#include "FlrwTestDataMixin.h"

#include <cppunit/extensions/HelperMacros.h>
#include <stdexcept>

class FlrwNatTest : public CppUnit::TestFixture, public FlrwTestDataMixin
{
    CPPUNIT_TEST_SUITE(FlrwNatTest);
    CPPUNIT_TEST_EXCEPTION(testMatterLessThanZeroThrows, std::domain_error);
    CPPUNIT_TEST_EXCEPTION(testVacuumLessThanZeroThrows, std::domain_error);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows21, std::domain_error);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows22, std::domain_error);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows23, std::domain_error);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows24, std::domain_error);
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

    /** Tests matter density < 0 */
    void testMatterLessThanZeroThrows();

    /** Tests vacuum energy density < 0 (Universe recollapses) */
    void testVacuumLessThanZeroThrows();

    /** Tests no Big Bang where om < 0.5 and b = 2 */
    void testNoBigBangThrows21();
    /** Tests no Big Bang where om < 0.5 and b < 2 */
    void testNoBigBangThrows22();
    /** Tests no Big Bang where om > 0.5 and b = 2 */
    void testNoBigBangThrows23();
    /** Tests no Big Bang where om > 0.5 and b = 2 */
    void testNoBigBangThrows24();

    void testLuminosityDistance();

    void testComovingDistance();

    void testComovingTransverseDistance();

    void testAngularDistance();

    void testAge();

    void testComovingVolume();
};


#endif // MILIA_FLRW_NAT_TEST_H
