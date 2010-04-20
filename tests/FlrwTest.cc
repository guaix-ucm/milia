/*
 * Copyright 2008-2010 Sergio Pascual
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

#include "FlrwTest.h"
#include "milia/metric.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FlrwTest);

void FlrwTest::setUp() {
}

void FlrwTest::tearDown() {
}

void FlrwTest::testHubbleZeroThrows() {
	const milia::metrics::flrw test00(0, 1, 1);
}

void FlrwTest::testHubbleLessThanZeroThrows() {
	const milia::metrics::flrw test00(-50, 1, 1);
}

void FlrwTest::testMatterLessThanZeroThrows() {
	const milia::metrics::flrw test00(50, -1, 1);
}

void FlrwTest::testVacuumLessThanZeroThrows() {
	const milia::metrics::flrw test00(50, 1, -1); // Recollapse
}

void FlrwTest::testRecollapse11Throws() {
	const milia::metrics::flrw test00(50, 1.5, 0.008665856); // Recollapse b = 2
}

void FlrwTest::testRecollapse12Throws() {
	const milia::metrics::flrw test00(50, 1.5, 0.007); // Recollapse b < 2
}

void FlrwTest::testNoBigBangThrows21() {
	const milia::metrics::flrw test00(50, 0.3, 1.713460403); // No Big Bang, b = 2 and om < 0.5
}

void FlrwTest::testNoBigBangThrows22() {
	const milia::metrics::flrw test00(50, 0.3, 2); // No Big Bang, b < 2 and om < 0.5
}

void FlrwTest::testNoBigBangThrows23() {
	const milia::metrics::flrw test00(50, 0.7, 2.254425343); // No Big Bang, b = 2 and om > 0.5
}

void FlrwTest::testNoBigBangThrows24() {
	const milia::metrics::flrw test00(50, 0.7, 3); // No Big Bang, b < 2 and om > 0.5
}

void FlrwTest::testLuminosityDistance() {
	/* test cases
	 * O_m == 0 and O_v == 0               case OM_OV_0
	 * O_m == 0 and O_v < 1                case OM
	 * O_m == 0 and O_v == 0               case OM_EDS
	 * 0 < O_m < 1 and O_v == 0            case OV_1
	 * O_m > 1 and O_v == 0                case OV_2
	 * O_m == 1 and O_v == 0               case OV_DS
	 * O_v != 0 and O_m + O_v == 1         case OM_OV_1
	 * O_m + O_v != 1 and b == 2           case A2_1
	 * O_m + O_v != 1 and (b > 2 or b < 0) case A1
	 * O_m + O_v != 1 and 0 < b < 2        case A2_2
	 */

	// Number of lum_models
	const int val = 7;
	for (int j = 0; j < val; ++j) {
		const milia::metrics::flrw test00(lum_model[j][0], lum_model[j][1],
				lum_model[j][2]);
		for (int i = 0; i < 5; ++i) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(lum_table[j][i][0],
					test00.dl(lum_table[j][i][1]), lum_table[j][i][2]);
		}
	}
}

void FlrwTest::testAngularDistance() {
	// Number of ang_models
	const int val = 1;
	for (int j = 0; j < val; ++j) {
		const milia::metrics::flrw test00(ang_model[j][0], ang_model[j][1],
				ang_model[j][2]);
		for (int i = 0; i < 5; ++i) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(ang_table[j][i][0],
					test00.da(ang_table[j][i][1]), ang_table[j][i][2]);
		}
	}
}

void FlrwTest::testComovingTransverseDistance() {
	// Number of cotran_models
	const int val = 1;
	for (int j = 0; j < val; ++j) {
		const milia::metrics::flrw test00(cotran_model[j][0], cotran_model[j][1],
				cotran_model[j][2]);
		for (int i = 0; i < 5; ++i) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(cotran_table[j][i][0],
					test00.dm(cotran_table[j][i][1]), cotran_table[j][i][2]);
		}
	}
}
void FlrwTest::testComovingDistance() {
	// Number of com_models
	const int val = 3;
	for (int j = 0; j < val; ++j) {
		const milia::metrics::flrw test00(com_model[j][0], com_model[j][1],
				com_model[j][2]);
		for (int i = 0; i < 5; ++i) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(com_table[j][i][0],
					test00.dc(com_table[j][i][1]), com_table[j][i][2]);
		}
	}
}

void FlrwTest::testAge() {
	/* test cases
	 * O_m == 0 and O_v == 0               case OM_OV_0
	 * O_m == 0 and O_v < 1                case OM
	 * O_m == 0 and O_v == 1               case OM_EDS
	 * 0 < O_m < 1 and O_v == 0            case OV_1
	 * O_m > 1 and O_v == 0                case OV_2
	 * O_m == 1 and O_v == 0               case OV_DS
	 * O_v != 0 and O_m + O_v == 1         case OM_OV_1
	 * O_m + O_v != 1 and b == 2           case A2_1
	 * O_m + O_v != 1 and (b > 2 or b < 0) case A1
	 * O_m + O_v != 1 and 0 < b < 2        case A2_2
	 */

	// Number of age_models
	const int val = 7;
	for (int j = 0; j < val; ++j) {
		const milia::metrics::flrw test00(age_model[j][0], age_model[j][1],
				age_model[j][2]);
		for (int i = 0; i < 5; ++i) {
			CPPUNIT_ASSERT_DOUBLES_EQUAL(age_table[j][i][0],
					test00.age(age_table[j][i][1]), age_table[j][i][2]);
		}
	}
}

void FlrwTest::testComovingVolume() {
	// Number of vol_models
	const int val = 4;
	for (int j = 0; j < val; ++j) {
		const milia::metrics::flrw test00(vol_model[j][0], vol_model[j][1],
				vol_model[j][2]);
		for (int i = 0; i < 5; ++i) {			
			CPPUNIT_ASSERT_DOUBLES_EQUAL(vol_table[j][i][0],
					test00.vol(vol_table[j][i][1]), vol_table[j][i][2]);
		}
	}
}
