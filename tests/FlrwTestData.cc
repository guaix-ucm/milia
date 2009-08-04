/*
 * Copyright 2008-2009 Sergio Pascual
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

#include "FlrwTestData.h"

   /* These tests are removed
     * {50.,0.,-1.} : recollapse
     * {50., 0.,1.5} : recollapse
     * {50., 1.5, 0.} : no Big Bang
     * A2_1 A2_2 recollapse
     * O_m + O_v != 1 and (b > 2 or b < 0) case A1
     * O_m + O_v != 1 and 0 < b < 2        case A2_2
     */

const double FlrwTestData::lum_model[][3] = {
    {50., 0.0, 0.0}, // OM_OV_0
    {70., 0.3, 0.7}, // OM_OV_1
    {50., 0.0, 0.5}, // OM
    {50., 0.5, 0.}, // OV_1
    //{50., 1.5, 0.}, // OV_2 makes the Universe recollapse
    {50., 1.0, 0.}, // OV_3
    {50., 0.3, 0.2}, // A1
    //{50., 1.5, 0.008665856}, // A2_1 recollapse
    //{50., 1.5, 0.007},    // A2_2
};

// {50., 0.3, 1.713460403} breaks
// {50., 0.7, 2.254425343} exception
// {50., 0.7, 3} exception
// {50., 0.3, 2.} exception

const double FlrwTestData::lum_table[][5][3] =
{
    {{0.0010005, 0.001, 1e-7}, // OM_OV_0
     {0.01005, 0.01, 1e-6},
     {0.105, 0.1, 1e-3},
     {1.5, 1, 1e-1},
     {60, 10, 1e-3}},
    {{4.28606, 0.001, 1e-5}, // OM_OV_1
     {43.1582, 0.01, 1e-4},
     {460.299, 0.1, 1e-3},
     {6607.65, 1, 1e-2},
     {103843, 10, 1}},
    {{0, 0.001, 1e-5}, // OM
     {0, 0.01, 1e-4},
     {0, 0.1, 1e-3},
     {0, 1, 1e-2},
     {0, 10, 1e-1}},
    {{5.99809, 0.001, 1e-5}, // OV_1
     {60.1827, 0.01, 1e-4},
     {621.523, 0.1, 1e-3},
     {7812.95, 1, 1e-2},
     {135543., 10, 1}},
/*    {{5.99659, 0.001, 1e-5}, // OV_2
     {60.0328, 0.01, 1e-4},
     {606.564, 0.1, 1e-3},
     {6445.82, 1, 1e-2},
     {71950.1, 10, 1e-1}},*/
    {{5.99734, 0.001, 1e-5}, // OV_3
     {60.1076, 0.01, 1e-4},
     {613.868, 0.1, 1e-3},
     {7024.56, 1, 1e-2},
     {92136.6, 10, 1e-1}},
    {{5.99899, 0.001, 1e-5}, // A1
     {60.2721, 0.01, 1e-4},
     {630.080, 0.1, 1e-3},
     {8470.21, 1, 1e-2},
     {168260., 10, 1}},
 /*   {{5.99661, 0.001, 1e-5}, // A2_1
     {60.0353, 0.01, 1e-4},
     {606.769, 0.1, 1e-3},
     {6449.59, 1, 1e-2},
     {71827.0, 10, 1}},*/
/*    {{5.99661, 0.001, 1e-5}, // A2_2
     {60.0348, 0.01, 1e-4},
     {606.730, 0.1, 1e-3},
     {6448.86, 1, 1e-2},
     {71850.7, 10, 1e-1}},*/
};

const double FlrwTestData::ang_model[][3] = {
    {50., 1.0, 0.},
    };

const double FlrwTestData::ang_table[][5][3] =
  {{{5.98536, 0.001, 1e-5},
    {58.9232, 0.01, 1e-4},
    {507.329, 0.1, 1e-3},
    {1756.14, 1, 1e-2},
    {761.459, 10, 1e-3}}};

const double FlrwTestData::cotran_model[][3] = {
    {50., 1.0, 0.},
    };

const double FlrwTestData::cotran_table[][5][3] =
  {{{5.99135, 0.001, 1e-5},
    {59.5124, 0.01, 1e-4},
    {558.062, 0.1, 1e-3},
    {3512.28, 1, 1e-2},
    {8376.05, 10, 1e-1}}};

const double FlrwTestData::com_model[][3] = {
    {50., 1.0, 0.},
    {50., 0.5, 0.6},
    {50., 0.5, 0.4},
    };

const double FlrwTestData::com_table[][5][3] =
  {
      {{5.99135, 0.001, 1e-5},
       {59.5124, 0.01, 1e-4},
       {558.062, 0.1, 1e-3},
       {3512.28, 1, 1e-2},
       {8376.05, 10, 1e-1}},
      {{5.99389, 0.001, 1e-5},
       {59.7634, 0.01, 1e-4},
       {580.027, 0.1, 1e-3},
       {4275.74, 1, 1e-2},
       {11212.1, 10, 1e-1}},
     {{5.99329, 0.001, 1e-5},
       {59.7041, 0.01, 1e-4},
       {574.702, 0.1, 1e-3},
       {4085.51, 1, 1e-2},
       {10714.2, 10, 1e-1}}
  };

const double FlrwTestData::age_model[][3] = {
    {50., 0.0, 0.0}, // OM_OV_0
    {70., 0.3, 0.7}, // OM_OV_1
    {50., 0.0, 0.5}, // OM
    {50., 0.5, 0.}, // OV_1
   // {50., 1.5, 0.}, // OV_2
    {50., 1.0, 0.}, // OV_3
    {50., 0.3, 0.2}, // A1
  //  {50., 1.5, 0.008665856}, // A2_1
  //  {50., 1.5, 0.007},    // A2_2
};

// {50., 0.3, 1.713460403} breaks
// {50., 0.7, 2.254425343} exception
// {50., 0.7, 3} exception
// {50., 0.3, 2.} exception

const double FlrwTestData::age_table[][5][3] =
{
 {{19.56, 0,1e-2}, // OM_OV_0
  {17.7818, 0.1, 1e-4},
  {9.78, 1, 1e-2},
  {1.77818, 10, 1e-5},
  {0.193663, 100, 1e-5}},
 {{13.4698, 0, 1e-4}, // OM_OV_1
  {12.1683, 0.1, 1e-4},
  {5.75287, 1, 1e-5},
  {0.465986, 10, 1e-6},
  {0.0167535, 100, 1e-7}},
 {{24.3806, 0, 1e-4}, // OM
  {22.5614, 0.1, 1e-4},
  {13.3113, 1, 1e-4},
  {2.51128, 10, 1e-5},
  {0.273877, 100, 1e-6}},
 {{14.7394, 0, 1e-4},  // OV_1
  {12.9823, 0.1, 1e-4},
  {5.74115, 1, 1e-5},
  {0.492328, 10, 1e-6},
  {0.0181145, 100, 1e-7}},
/* {{11.9562, 0, 1e-4},  // OV_2
  {10.2382, 0.1, 1e-4},
  {3.97141, 1, 1e-5},
  {0.294536, 10, 1e-6},
  {0.0104998, 100, 1e-7}},*/
 {{13.04, 0, 1e-2},  // OV_3
  {11.3029, 0.1, 1e-4},
  {4.61034, 1, 1e-5},
  {0.357428, 10, 1e-6},
  {0.0128468, 100, 1e-7}},
 {{21.4336, 0, 1e-4}, // A1
  {18.8297, 0.1, 1e-4},
  {8.29440, 1, 1e-5},
  {0.736001, 10, 1e-6},
  {0.0280273, 100, 1e-7}},
/* {{313.408, 0, 1e-3}, // A2_1
  {295.759, 0.1, 1e-3},
  {209.899, 1, 1e-3},
  {85.938, 10, 1e-3},
  {28.1435, 100, 1e-4}},
 {{-41.3425, 0, 1e-4}, // A2_2
  {-36.0076, 0.1, 1e-4},
  {-15.1328, 1, 1e-4},
  {-1.22642, 10, 1e-5},
  {-0.0445536, 100, 1e-7}}*/
};

const double FlrwTestData::vol_model[][3] = {
    {50., 0.5, 0.},
    //{50., 1.5, 0.},
    {50., 1.0, 0.},
};

// {50., 0.3, 1.713460403} breaks
// {50., 0.7, 2.254425343} exception
// {50., 0.7, 3} exception
// {50., 0.3, 2.} exception

const double FlrwTestData::vol_table[][5][3] =
{
  {
	  {89304527.4589468, 1e-3, 1e-1},
	  {8.94843e+08, 1e-2, 1e3},
	  {9.09771e+09, 1e-1, 1e4},
	  {9.27905e+10, 1 , 1e5},
	  {5.35838e+11, 10, 1e6},
  },
 /* {
	  {89434230.119928,1e-3,1e-1},
	  {9.07466e+08,1e-2,1e3},
	  {1.00685e+10,1e-1,1e5},
	  {101367659891.337,1,1e5},
	  {2.71387e+11,10,1e6}
  },*/
  {
		  {71.6889,1e-3,1e-4},
		  {70259,1e-2,1},
		  {5.79329e+07,1e-1,1e2},
		  {1.44426e+10,1,1e5},
		  {1.95883e+11,10,1e6},
  },

};
