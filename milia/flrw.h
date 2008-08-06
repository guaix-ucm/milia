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

#ifndef _MILIA_FLRW_H_
#define _MILIA_FLRW_H_

namespace milia {

namespace metrics {

/**
 * The Friedmann-Lemaître-Robertson-Walker metric
 */
class flrw {
public:
	flrw(double hubble_constant, double matter_density, double lambda_density);

	double distance_luminosity(double redshift) const;
	double distance_comoving(double redshift) const;
	double distance_angular(double redshift) const;
	double time_lookback(double redshift) const;

private:
	//!Hubble radius
	static const double ms_hubble_radius;
	//!Hubble time
	static const double ms_hubble_time;

	//!Hubble constant
	double m_hu;
	//!Matter density
	double m_om;
	//!cosmologyologycal constant density
	double m_ol;

};

} // namespace metrics

} // namespace milia


#endif /* _MILIA_METRIC_H_ */
