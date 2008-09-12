/*
 * Copyright 2008 Sergio Pascual
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

#include "FlrwTest.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FlrwTest);

void FlrwTest::setUp()
{}

void FlrwTest::tearDown()
{}

void FlrwTest::testHubbleZeroThrows()
{
    const milia::metric test00(0, 1, 1);
}

void FlrwTest::testHubbleLessThanZeroThrows()
{
    const milia::metric test00(-50, 1, 1);
}

void FlrwTest::testMatterLessThanZeroThrows()
{
    const milia::metric test00(50, -1, 1);
}

void FlrwTest::testVacuumLessThanZeroThrows()
{
    const milia::metric test00(50, 1, -1); // Recollapse
}

void FlrwTest::testRecollapse11Throws()
{
    const milia::metric test00(50, 1.5, 0.008665856); // Recollapse b = 2
}

void FlrwTest::testRecollapse12Throws()
{
    const milia::metric test00(50, 1.5, 0.007); // Recollapse b < 2
}

void FlrwTest::testNoBigBangThrows21()
{
    const milia::metric test00(50, 0.3, 1.713460403); // No Big Bang, b = 2 and om < 0.5
}

void FlrwTest::testNoBigBangThrows22()
{
    const milia::metric test00(50, 0.3, 2); // No Big Bang, b < 2 and om < 0.5
}

void FlrwTest::testNoBigBangThrows23()
{
    const milia::metric test00(50, 0.7, 2.254425343); // No Big Bang, b = 2 and om > 0.5
}

void FlrwTest::testNoBigBangThrows24()
{
    const milia::metric test00(50, 0.7, 3); // No Big Bang, b < 2 and om > 0.5
}

void FlrwTest::testLuminosityDistance()
{
    const milia::metric test00(50, 1, 0);
    const int val = 5;
    const double tabulated_values[val][3] = {{5.99734, 0.001, 1e-5},                                            
                                            {60.1076, 0.01, 1e-4},
                                            {613.868, 0.1, 1e-3},
                                            {7024.56, 1, 1e-2},
                                            {92136.6, 10, 1e-1}};

    for(int i = 0; i < val; ++i)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(tabulated_values[i][0],
                                     test00.distance_luminosity(tabulated_values[i][1]),
                                     tabulated_values[i][2]);    
}




