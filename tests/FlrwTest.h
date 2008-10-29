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

#ifndef MILIA_FLRW_TEST_H
#define MILIA_FLRW_TEST_H

#include "milia/exception.h"

#include <cppunit/extensions/HelperMacros.h>

class FlrwTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(FlrwTest);
 /*   CPPUNIT_TEST_EXCEPTION(testHubbleZeroThrows, milia::exception);
    CPPUNIT_TEST_EXCEPTION(testHubbleLessThanZeroThrows, milia::exception);
    CPPUNIT_TEST_EXCEPTION(testMatterLessThanZeroThrows, milia::exception);
    CPPUNIT_TEST_EXCEPTION(testVacuumLessThanZeroThrows, milia::recollapse);
    CPPUNIT_TEST_EXCEPTION(testRecollapse11Throws, milia::recollapse);
    CPPUNIT_TEST_EXCEPTION(testRecollapse12Throws, milia::recollapse);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows21, milia::no_big_bang);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows22, milia::no_big_bang);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows23, milia::no_big_bang);
    CPPUNIT_TEST_EXCEPTION(testNoBigBangThrows24, milia::no_big_bang);*/
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

    /** Tests matter density < 0 */
    void testMatterLessThanZeroThrows();

    /** Tests vacuum energy density < 0 (Universe recollapses) */
    void testVacuumLessThanZeroThrows();
    /** Tests recollapse where ov > 0 and b = 2 */
    void testRecollapse11Throws();
    /** Tests recollapse where ov > 0 and b < 2 */
    void testRecollapse12Throws();

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

private:
    static const double lum_model[][3];
    static const double lum_table[][5][3];
    static const double ang_model[][3];
    static const double ang_table[][5][3];
    static const double cotran_model[][3];
    static const double cotran_table[][5][3];
    static const double com_model[][3];
    static const double com_table[][5][3];
    static const double age_model[][3];
    static const double age_table[][5][3];
    static const double vol_model[][3];
    static const double vol_table[][5][3];
};


#endif // MILIA_FLRW_TEST_H
