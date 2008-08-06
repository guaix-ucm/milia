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

#include "flrw.h"
#include "exception.h"

namespace milia {

namespace metrics {

flrw::flrw(double hubble, double matter, double lambda) :
	m_hu(hubble), m_om(matter), m_ol(lambda)

{
	// om < 0
	if (m_om < 0)
		throw milia::exception();
	//ol < 0
	if (m_ol < 0)
		throw milia::exception();
}

double flrw::distance_luminosity(double redshift) const {
	return 0;
}

double flrw::distance_comoving(double redshift) const {
	return 0;
}

double flrw::distance_angular(double redshift) const {
	return 0;
}

double flrw::time_lookback(double redshift) const {
	return 0;
}

} // namespace metrics

} // namespace milia