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
#include "milia/flrw_nat.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FlrwTestNew);

void FlrwTestNew::setUp() {
}

void FlrwTestNew::tearDown() {
}

void FlrwTestNew::testMatterLessThanZeroThrows() {
  const milia::flrw_nat_new test00(-1, 1);
} 

void FlrwTestNew::testVacuumLessThanZeroThrows() {
   const milia::flrw_nat_new test00(1, -1); // Recollapse
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
        //const int val = 7;
        const int val = 2;
        for (int j = 1; j < val; ++j) {
          const double h = lum_model[j][0];
          const double m = lum_model[j][1];
          const double v = lum_model[j][2];
          const milia::flrw_nat_new test00(m, v);
          for (int i = 0; i < 5; ++i) {
            const double dis = lum_table[j][i][0] * h / 299792.458;
            const double z = lum_table[j][i][1];
            const double tol = lum_table[j][i][2];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dis, test00.dl(z), tol);
          }
        }
  }

  void FlrwTestNew::testComovingDistance() {
        const int val = 1;
        for (int j = 0; j < val; ++j) {
          const double h = com_model[j][0];
          const double m = com_model[j][1];
          const double v = com_model[j][2];
          const milia::flrw_nat_new test00(m, v);
          for (int i = 0; i < 5; ++i) {
            const double dis = com_table[j][i][0] * h / 299792.458;
            const double z = com_table[j][i][1];
            const double tol = com_table[j][i][2];
            CPPUNIT_ASSERT_DOUBLES_EQUAL(dis, test00.dc(z), tol);
          }
        }
   }

  void FlrwTestNew::testComovingTransverseDistance() {
        const int val = 0;
        for (int j = 0; j < val; ++j) {
          const milia::flrw_nat_new test00(cotran_model[j][0], cotran_model[j][1]);
            for (int i = 0; i < 5; ++i) {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(cotran_table[j][i][0],
                   test00.dm(cotran_table[j][i][1]), cotran_table[j][i][2]);
            }
        }
   }

  void FlrwTestNew::testAngularDistance() {
        const int val = 1;
        for (int j = 0; j < val; ++j) {
          const milia::flrw_nat_new test00(ang_model[j][0], ang_model[j][1]);
            for (int i = 0; i < 5; ++i) {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(cotran_table[j][i][0],
                   test00.da(ang_table[j][i][1]), ang_table[j][i][2]);
            }
        }
   }