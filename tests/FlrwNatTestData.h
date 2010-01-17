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

#ifndef MILIA_FLRW_NAT_TEST_DATA_H
#define MILIA_FLRW_NAT_TEST_DATA_H

struct FlrwNatTestData
{
    static const double lum_model[][2];
    static const double lum_table[][5][3];
    static const double ang_model[][2];
    static const double ang_table[][5][3];
    static const double cotran_model[][2];
    static const double cotran_table[][5][3];
    static const double com_model[][2];
    static const double com_table[][5][3];
    static const double age_model[][2];
    static const double age_table[][5][3];
    static const double vol_model[][2];
    static const double vol_table[][5][3];
};


#endif // MILIA_FLRW_NAT_TEST_DATA_H
