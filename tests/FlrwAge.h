/*
 * Copyright 2009 Sergio Pascual
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

#ifndef MILIA_FLRW_AGE_H
#define MILIA_FLRW_AGE_H

#include "milia/exception.h"

#include <cppunit/extensions/HelperMacros.h>

class FlrwAge : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(FlrwAge);
    CPPUNIT_TEST(testAge);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp();

    void tearDown();

    /** Checks the age is OK */
    void testAge();
};


#endif // MILIA_FLRW_AGE_H
