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
#include <cmath>
#include <gsl/gsl_math.h>

namespace milia
{

    namespace metrics
    {

        const double flrw::ms_hubble_radius = 2.99792e5; // Hubble Radius in Mpc
        const double flrw::ms_hubble_time = 9.78e2; // Hubble time in Gyr


        flrw::flrw(double hubble, double matter, double vacuum) :
                m_hu(hubble), m_om(matter), m_ov(vacuum), m_ok(1.0 - matter - vacuum),
                m_kap(GSL_SIGN(m_ok))

        {
            // hubble <= 0
            if (m_hu <= 0)
                throw milia::exception("Hubble constant <= 0 not allowed");
            // om < 0
            if (m_om < 0)
                throw milia::exception("Matter density < 0 not allowed");
            //ol < 0
            if (m_ov < 0)
                throw milia::recollapse("The Universe recollapses"); // Recollapse

            m_b = (-13.5 * gsl_pow_2(m_om) * m_ov / (gsl_pow_3(m_ok)));

            if(m_om >= 1 and m_b <= 2)
                throw milia::recollapse("The Universe recollapses"); // Recollapse

            if(m_ov >= 1 and m_b <= 2)
                throw milia::no_big_bang("No Big Bang"); // No Big bang with this parameters
            
            m_case = calculation_case(m_om, m_ov, m_ok, m_b);

            m_r_h = ms_hubble_radius / m_hu;
            m_t_h = ms_hubble_time / m_hu;
            m_universe_age = age();
        }

        flrw::CALC_CASES flrw::calculation_case(double matter, double lambda,
                                                double curvature, double b)
        {
            const bool l3 = (matter <= 0);
            const bool l4 = (lambda <= 0);
            // om=ol=0
            if (l3 && l4)
                return OM_OV_0;
            // ol=0 0<om<1
            if (l4 && matter < 1)
                return OV_1;
            // ol=0 om>1
            if (l4 && matter> 1)
                return OV_2;
            // ol=0 om==1
            if (l4)
                return OV_3;
            // om=0
            if (l3)
                return OM;
            // om+ol==1
            if (curvature <= 0)
                return OM_OV_1;
            // om+ol!=1
            // a1 b<0 || b>2
            if ((b <= 0) ||(b> 2))
                return A1;
            // a2 b=2
            if (b == 2)
                return A2_1;
            // a2 0<b<2
            if ((b> 0) && (b < 2))
                return A2_2;
            // Fallback
            return NO_CASE;

        }

        double flrw::helper(double z) const
        {
            return sqrt(m_om*gsl_pow_3(1+z)+m_ok*gsl_pow_2(1+z)+m_ov);
            // No se cual de estas dos es mas rapida
            // return sqrt(gsl_pow_2(1.+z)*(1.+m_om*z)-z*m_ol*(2.+z))));
        }

        double flrw::distance_comoving(double z) const
        {
            const double dm = distance_comoving_transverse(z);
            if (m_ok == 0)
                return dm;
            else {
                const double sqok = sqrt(std::abs(m_ok));                
                if (m_ok > 0)
                    return m_r_h / sqok * asinh(dm * sqok / m_r_h);
                else
                    return m_r_h / sqok * asin(dm * sqok / m_r_h);
            }
        }

        double flrw::distance_comoving_transverse(double z) const
        {
            return distance_luminosity(z) / (1. + z);
        }

        double flrw::distance_angular(double z) const
        {
            return distance_luminosity(z) / gsl_pow_2(1. + z);
        }

        double flrw::lookback_time(double z) const
        {
            return m_universe_age - age(z);
        }

    } // namespace metrics

}
// namespace milia
