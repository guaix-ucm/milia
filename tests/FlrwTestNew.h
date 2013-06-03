/*
 * Copyright 2008-2013 Sergio Pascual
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

#ifndef MILIA_FLRW_TEST_NEW_H
#define MILIA_FLRW_TEST_NEW_H

#include "FlrwTestData.h"
#include <stdexcept>
#include <cppunit/extensions/HelperMacros.h>

class FlrwTestNew : public CppUnit::TestFixture, public FlrwTestData
{
    CPPUNIT_TEST_SUITE(FlrwTestNew);
    CPPUNIT_TEST_EXCEPTION(testMatterLessThanZeroThrows, std::domain_error);
    CPPUNIT_TEST_EXCEPTION(testVacuumLessThanZeroThrows, std::domain_error);
    CPPUNIT_TEST(testLuminosityDistance);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();

    void tearDown();

    void testMatterLessThanZeroThrows();
    void testVacuumLessThanZeroThrows();
    void testLuminosityDistance();
};


#endif // MILIA_FLRW_TEST_NEW_H
