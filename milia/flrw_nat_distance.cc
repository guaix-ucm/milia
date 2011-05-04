/*
 * Copyright 2008-2011 Sergio Pascual
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "flrw_nat.h"
#include "flrw_prec.h"
#include "exception.h"

#include <cmath>
#include <cstdlib>

// Provides M_SQRT3
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_ellint.h>

namespace
{
  static const double M_4THRT3 = sqrt(M_SQRT3);
  static const double FLAT_K = sqrt(0.5 + 0.25 * M_SQRT3);
}

namespace milia
{
  namespace metrics
  {
    using std::abs;

    // Luminosity distance
    double flrw_nat::dl(double z) const
    {
      switch (m_case)
      {
        case OM_OV_0:
          return 0.5 * z * (z + 2);
          break;
        case OV_1:
        case OV_2:
        case OV_EDS:
          return 2. * ((2. - m_om * (1. - z) - (2. - m_om) * sqrt(1
              + m_om * z))) / (m_om * m_om);
          break;
        case OM:
          return ((1 + z) / m_ov) * (1 + z - sqrt(m_ov + (1 - m_ov)
              * gsl_pow_2(1 + z)));
        case OM_DS:
          return z * (1 + z);
          break;
        case A1:
          //om+ol != 1, crit < 0 || crit > 2
        {
          const double v = pow(m_kap * (m_crit - 1.) + sqrt(m_crit * (m_crit
              - 2.)), 1. / 3.);
          const double y = (-1. + m_kap * (v + 1. / v)) / 3.;
          const double A = sqrt(y * (3. * y + 2.));
          const double g = 1. / sqrt(A);
          const double k = sqrt(0.5 + 0.25 * gsl_pow_2(g) * (v + 1. / v));
          const double sup = m_om / abs(m_ok);
          const double phi = acos(((1. + z) * sup + m_kap * y - A) / ((1. + z)
              * sup + m_kap * y + A));
          const double phi0 = acos((sup + m_kap * y - A)
              / (sup + m_kap * y + A));
          return (1 + z) / m_sqok * sinc(m_kap, 1.0, g
              * (gsl_sf_ellint_F(phi0, k, ELLIP_PREC) - gsl_sf_ellint_F(phi, k,
                  ELLIP_PREC)));
        }
        case A2_1:
          // crit = 2, lower region
        case A2_2:
          // 0 < crit < 2, lower region
        {
          const double arg0 = acos(1. - m_crit) / 3.;
          const double arg1 = m_om / m_sqok;
          const double y1 = (-1. + cos(arg0) + M_SQRT3 * sin(arg0)) / 3.;
          const double y2 = (-1. - 2. * cos(arg0)) / 3.;
          const double y3 = (-1. + cos(arg0) - M_SQRT3 * sin(arg0)) / 3.;
          const double g = 2. / sqrt(y1 - y2);
          const double k = sqrt((y1 - y3) / (y1 - y2));
          const double phi = asin(sqrt((y1 - y2) / ((1 + z) * arg1 + y1)));
          const double phi0 = asin(sqrt((y1 - y2) / (arg1 + y1)));
          return (1. + z) / m_sqok * sin(g * (gsl_sf_ellint_F(phi0, k,
              ELLIP_PREC) - gsl_sf_ellint_F(phi, k, ELLIP_PREC)));
        }
        case OM_OV_1:
        {
          // om + ov = 1
          const double arg0 = pow((1 / m_om - 1), 1. / 3);
          const double down = 1 + (1 + M_SQRT3) * arg0;
          const double up = 1 + (1 - M_SQRT3) * arg0;
          const double phi = acos((z + up) / (z + down));
          const double phi0 = acos(up / down);
          return (1 + z) / (M_4THRT3 * sqrt(m_om) * sqrt(arg0))
              * (gsl_sf_ellint_F(phi0, FLAT_K, ELLIP_PREC) - gsl_sf_ellint_F(
                  phi, FLAT_K, ELLIP_PREC));
        }
        default:
          break;
      }
      return -1;
    }
  } // namespace metrics
} //namespace milia
