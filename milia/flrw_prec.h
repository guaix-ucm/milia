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

#ifndef MILIA_FLRW_PREC_H
#define MILIA_FLRW_PREC_H

#define FLRW_EQ_TOL 1.e-14

#include <gsl/gsl_mode.h>
/* Precision of the elliptical funcions in gsl */
/* Single is about 10^-7 */
#define ELLIP_PREC GSL_PREC_SINGLE

#endif /* MILIA_FLRW_PREC_H */
