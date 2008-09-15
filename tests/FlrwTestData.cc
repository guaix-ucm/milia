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

  /* test cases for luminosity distance
   * O_m == 0 and O_v == 0               case OM_OV_0
   * O_m == 0 and O_v != 0               case OM
   * 0 < O_m < 1 and O_v == 0            case OV_1
   * O_m > 1 and O_v == 0                case OV_2
   * O_m == 1 and O_v == 0               case OV_3
   * O_v != 0 and O_m + O_v == 1         case OM_OV_1
   * O_m + O_v != 1 and b == 2           case A2_1
   * O_m + O_v != 1 and (b > 2 or b < 0) case A1
   * O_m + O_v != 1 and 0 < b < 2        case A2_2
   */
              
     
 
/*
    const double lum_model[][3] = {
        {50., 0.,0.}, // OM_OV_0 
        {50.,0.,0.5}, {50.,0.,1.}, //OM 
        {50., 0.5, 0.}, {50., 1.0, 0.}, // OV_1 OV_2 OV_3
        {70, 0.3, 0.7}, {50., 0., 1.} // OM_OV_1        
        };
    */

   /* These tests are removed
     * {50.,0.,-1.} : recollapse
     * {50., 0.,1.5} : recollapse
     * {50., 1.5, 0.} : no Big Bang
     */


const double FlrwTest::lum_model[][3] = {             
    {50., 1.0, 0.},  
    };

const double FlrwTest::lum_table[][5][3] = 
  {{{5.99734, 0.001, 1e-5},                                            
    {60.1076, 0.01, 1e-4},
    {613.868, 0.1, 1e-3},
    {7024.56, 1, 1e-2},
    {92136.6, 10, 1e-1}}};
