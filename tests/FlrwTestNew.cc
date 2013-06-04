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

#include "FlrwTestNew.h"
#include "milia/flrw.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FlrwTestNew);

using milia::rei::flrw;

void FlrwTestNew::setUp() {
}

void FlrwTestNew::tearDown() {
}

void FlrwTestNew::testHubbleLessThanZeroThrows() {
  const flrw test00(-50, 1, 1);
}

void FlrwTestNew::testMatterLessThanZeroThrows() {
  const flrw test00(50, -1, 1);
} 

void FlrwTestNew::testVacuumLessThanZeroThrows() {
   const flrw test00(50, 1, -1); // Recollapse
}

void FlrwTestNew::testLuminosityDistance() {
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
          const double h = lum_model[j][0];
          const double m = lum_model[j][1];
          const double v = lum_model[j][2];
          const flrw test00(h, m, v);
          for (int i = 0; i < 5; ++i) {
            const double dis = lum_table[j][i][0];
            const double z = lum_table[j][i][1];
            const double tol = lum_table[j][i][2];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dis, test00.dl(z), tol);
          }
        }
  }

  void FlrwTestNew::testComovingDistance() {
        const int val = 3;
        for (int j = 0; j < val; ++j) {
          const double h = com_model[j][0];
          const double m = com_model[j][1];
          const double v = com_model[j][2];
          const flrw test00(h, m, v);
          for (int i = 0; i < 5; ++i) {
            const double dis = com_table[j][i][0];
            const double z = com_table[j][i][1];
            const double tol = com_table[j][i][2];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dis, test00.dc(z), tol);
          }
        }
   }

  void FlrwTestNew::testComovingTransverseDistance() {
        const int val = 1;
        for (int j = 0; j < val; ++j) {
          const double h = cotran_model[j][0];
          const double m = cotran_model[j][1];
          const double v = cotran_model[j][2];
          const flrw test00(h, m, v);
          for (int i = 0; i < 5; ++i) {
            const double dis = cotran_table[j][i][0];
            const double z = cotran_table[j][i][1];
            const double tol = cotran_table[j][i][2];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dis, test00.dm(z), tol);
            }
        }
   }

  void FlrwTestNew::testAngularDistance() {
        const int val = 1;
        for (int j = 0; j < val; ++j) {
          const double h = ang_model[j][0];
          const double m = ang_model[j][1];
          const double v = ang_model[j][2];
          const flrw test00(h, m, v);
          for (int i = 0; i < 5; ++i) {
            const double dis = ang_table[j][i][0];
            const double z = ang_table[j][i][1];
            const double tol = ang_table[j][i][2];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dis, test00.da(z), tol);
            }
        }
   }

void FlrwTestNew::testComovingVolume() {
        // Number of vol_models
        const int val = 4;
        for (int j = 0; j < val; ++j) {
                const flrw test00(vol_model[j][0], vol_model[j][1],
                                vol_model[j][2]);
                for (int i = 0; i < 5; ++i) {                   
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(vol_table[j][i][0],
                                        test00.vol(vol_table[j][i][1]), vol_table[j][i][2]);
                }
        }
}

void FlrwTestNew::testAge() {
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
                const flrw test00(age_model[j][0], age_model[j][1],
                                age_model[j][2]);
                for (int i = 0; i < 5; ++i) {
                        CPPUNIT_ASSERT_DOUBLES_EQUAL(age_table[j][i][0],
                                        test00.age(age_table[j][i][1]), age_table[j][i][2]);
                }
        }
}
