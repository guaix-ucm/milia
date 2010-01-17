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

// $Id$

#include <cmath>

#include "FlrwTest.h"
#include "milia/flrw.h"
#include "milia/flrw_nat.h"

using milia::metrics::flrw;
using milia::metrics::flrw_nat;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FlrwTest);


void FlrwTest::setUp() {
}

void FlrwTest::tearDown() {
}

void FlrwTest::testHubbleZeroThrows() {
	const milia::metrics::flrw phys(0, 1, 1);
}

void FlrwTest::testHubbleLessThanZeroThrows() {
	const milia::metrics::flrw phys(-50, 1, 1);
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
		const milia::metrics::flrw phys(ms_hubble, lum_model[j][0],
				lum_model[j][1]);
		const milia::metrics::flrw_nat nat(lum_model[j][0],
		        lum_model[j][1]);

		const double scale = phys.hubble_radius();
		const double rel = std::abs(scale * nat.dl(ms_z) / phys.dl(ms_z) - 1);
		CPPUNIT_ASSERT(rel < ms_rel_tol);
	}
}

void FlrwTest::testAngularDistance() {
	// Number of ang_models
	const int val = 1;
	for (int j = 0; j < val; ++j) {
		const flrw phys(ms_hubble, ang_model[j][0], ang_model[j][1]);
		const flrw_nat nat(ang_model[j][0], ang_model[j][1]);
		const double scale = phys.hubble_radius();
		const double rel = std::abs(scale * nat.da(ms_z) / phys.da(ms_z) - 1);
		CPPUNIT_ASSERT(rel < ms_rel_tol);
	}
}

void FlrwTest::testComovingTransverseDistance() {
	// Number of cotran_models
	const int val = 1;
	for (int j = 0; j < val; ++j) {
		const flrw phys(ms_hubble, cotran_model[j][0], cotran_model[j][1]);
		const flrw_nat nat(cotran_model[j][0], cotran_model[j][1]);
		const double scale = phys.hubble_radius();
		const double rel = std::abs(scale * nat.dm(ms_z) / phys.dm(ms_z) - 1);
		CPPUNIT_ASSERT(rel < ms_rel_tol);
	}
}
void FlrwTest::testComovingDistance() {
	// Number of com_models
	const int val = 3;
	for (int j = 0; j < val; ++j) {
		const flrw phys(ms_hubble, com_model[j][0], com_model[j][1]);
		const flrw_nat nat(com_model[j][0], com_model[j][1]);
		const double scale = phys.hubble_radius();
    const double rel = std::abs(scale * nat.dc(ms_z) / phys.dc(ms_z) - 1);
    CPPUNIT_ASSERT(rel < ms_rel_tol);
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
	const int val = 10;
	for (int j = 0; j < val; ++j) {
		const flrw phys(ms_hubble, age_model[j][0], age_model[j][1]);
		const flrw_nat nat(age_model[j][0], age_model[j][1]);
		const double scale = phys.hubble_time();
		const double rel = std::abs(scale * nat.age(ms_z) / phys.age(ms_z) - 1);
		CPPUNIT_ASSERT(rel < ms_rel_tol);
	}
}

void FlrwTest::testComovingVolume() {
	// Number of vol_models
	const int val = 4;
	for (int j = 0; j < val; ++j) {
		const flrw phys(ms_hubble, vol_model[j][0], vol_model[j][1]);
		const flrw_nat nat(vol_model[j][0], vol_model[j][1]);
		const double scale = phys.hubble_radius();
		const double vol = scale * scale * scale;
		const double rel = std::abs(vol * nat.vol(ms_z) / phys.vol(ms_z) - 1);
		CPPUNIT_ASSERT(rel < ms_rel_tol);
	}
}
