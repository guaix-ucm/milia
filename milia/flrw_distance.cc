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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <cstdlib>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp>

#include "flrw.h"

using std::sqrt;
using std::abs;
using std::sin;
using std::sinh;
using std::acos;
using std::atan;
using std::log;

using boost::math::ellint_1;
using boost::math::cbrt;
using boost::math::asinh;
using boost::math::atanh;
using boost::math::pow;

namespace
{
  const double M_SQRT3 = sqrt(3);
  const double M_4THRT3 = sqrt(M_SQRT3);
}

namespace milia
{
  namespace metrics
  {
    // Luminosity distance
    double flrw::dl(double z) const
    {
      switch (m_case)
      {
        case OM_OV_0:
          return m_r_h * 0.5 * z * (z + 2);
        case OV_1:
        case OV_2:
        case OV_EDS:
          return m_r_h * 2 * ((2 - m_om * (1 - z) - (2 - m_om) * sqrt(1 + m_om
              * z))) / pow<2> (m_om);
        case OM:
        case OM_DS:
        {
          // TODO: Implement this
          return 0;
        }
        case A1:
          //om+ol != 1 b < 0 || b > 2
        {
          const double v = cbrt(m_kap * (m_b - 1) + sqrt(m_b * (m_b - 2)));
          const double y = (-1 + m_kap * (v + 1. / v)) / 3.;
          const double A = sqrt(y * (3 * y + 2));
          const double g = 1. / sqrt(A);
          const double k = sqrt(0.5 + 0.25 * pow<2> (g) * (v + 1. / v));
          const double sup = m_om / abs(m_ok);
          const double phi = acos(((1 + z) * sup + m_kap * y - A) / ((1 + z)
              * sup + m_kap * y + A));
          const double phi0 = acos((sup + m_kap * y - A)
              / (sup + m_kap * y + A));
          return m_r_h * (1 + z) / m_sqok * sinc(m_kap, 1.0, g * (ellint_1(k,
              phi0) - ellint_1(k, phi)));
        }
        case A2_1: // b=2
        case A2_2: // 0 < b < 2
        {
          const double arg0 = acos(1 - m_b) / 3.;
          const double arg1 = m_om / abs(m_ok);
          const double y1 = (-1. + cos(arg0) + M_SQRT3 * sin(arg0)) / 3.;
          const double y2 = (-1. - 2. * cos(arg0)) / 3.;
          const double y3 = (-1. + cos(arg0) - M_SQRT3 * sin(arg0)) / 3.;
          const double g = 2. / sqrt(y1 - y2);
          const double k = sqrt((y1 - y3) / (y1 - y2));
          const double phi = asin(sqrt((y1 - y2) / ((1 + z) * arg1 + y1)));
          const double phi0 = asin(sqrt((y1 - y2) / (arg1 + y1)));
          return m_r_h * (1. + z) / m_sqok * sin(g * (ellint_1(k, phi0)
              - ellint_1(k, phi)));
        }
        case OM_OV_1:
        {
          // om + ol = 1
          const double k = sqrt(0.5 + 0.25 * M_SQRT3);
          const double c1 = M_4THRT3;
          const double arg0 = cbrt((1 / m_om - 1));
          const double down = 1 + (1 + M_SQRT3) * arg0;
          const double up = 1 + (1 - M_SQRT3) * arg0;
          const double phi = acos((z + up) / (z + down));
          const double phi0 = acos(up / down);
          return m_r_h * (1 + z) / (c1 * sqrt(m_om) * sqrt(arg0)) * (ellint_1(
              k, phi0) - ellint_1(k, phi));
        }
      }
      return -1;
    }
  } // namespace metrics
} // namespace milia
