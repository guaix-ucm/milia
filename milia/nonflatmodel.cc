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

#include "nonflatmodel.h"
#include "util.h"

#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/cbrt.hpp>

#include "flrw_nat_impl.h"
#include "util.h"

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
  namespace impl {

        double flrw_nat_nonflat::dc(double z) const {
          return asinc(m_kap, m_sqok, dm(z));
        }

        double flrw_nat_nonflat::vol(double z) const {
          const double lm = dm(z);
          return (lm * sqrt(1 + m_ok * pow<2> (lm)) - asinc(m_kap, m_sqok, lm)) / (2 * m_ok);
        }

        double flrw_nat_OM_OV_0::dl(double z) const 
        {
          return 0.5 * z * (z + 2);
        }

        double flrw_nat_OM_OV_0::age(double z) const
        {
          return 1.0 / (1 + z);
        }


      double flrw_nat_OV::dl(double z) const
      {
        return 2 * ((2 - m_om * (1 - z) - (2 - m_om) * sqrt(1 + m_om * z)))
                      / pow<2> (m_om);
      }

      double flrw_nat_OV_1::age(double z) const 
      {
        const double pre0 = 1 - m_om;
        const double prez = sqrt(1 + m_om * z);
        return (prez / (1 + z) - m_om / sqrt(pre0) * atanh(sqrt(pre0) / prez)) / pre0;
      }

      double flrw_nat_OV_2::age(double z) const 
      {
        const double pre0 = 1 - m_om;
        const double prez = sqrt(1 + m_om * z);
        return (prez / (1 + z) - m_om / sqrt(-pre0) * atan(sqrt(-pre0) / prez)) / pre0;
      }
   
   double flrw_nat_OM::dl(double z) const 
   {
     return ((1 + z) / m_ov) * (1 + z - sqrt(m_ov + (1 - m_ov)
                 * pow<2> (1 + z)));
   }

   double flrw_nat_OM::age(double z) const 
   {
     return asinh(1 / ((1 + z) * sqrt(1 / m_ov - 1))) / sqrt(m_ov);
   }

   double flrw_nat_A1::dl(double z) const
   {
      const double v = cbrt(m_kap * (m_crit - 1) + sqrt(m_crit * (m_crit - 2)));
      const double y = (-1 + m_kap * (v + 1. / v)) / 3.;
      const double A = sqrt(y * (3 * y + 2));
      const double g = 1. / sqrt(A);
      const double k = sqrt(0.5 + 0.25 * pow<2> (g) * (v + 1. / v));
      const double sup = m_om / abs(m_ok);
      const double phi = acos(((1 + z) * sup + m_kap * y - A) / ((1 + z)
          * sup + m_kap * y + A));
      const double phi0 = acos((sup + m_kap * y - A)
          / (sup + m_kap * y + A));
      return (1 + z) / m_sqok * sinc(m_kap, 1.0, g * (ellint_1(k,
          phi0) - ellint_1(k, phi)));
   }

   double flrw_nat_A1::age(double z) const
   {
     return 0;
   }

   double flrw_nat_A2::dl(double z) const
   {
    const double arg0 = acos(1 - m_crit) / 3.;
    const double arg1 = m_om / abs(m_ok);
    const double y1 = (-1. + cos(arg0) + M_SQRT3 * sin(arg0)) / 3.;
    const double y2 = (-1. - 2. * cos(arg0)) / 3.;
    const double y3 = (-1. + cos(arg0) - M_SQRT3 * sin(arg0)) / 3.;
    const double g = 2. / sqrt(y1 - y2);
    const double k = sqrt((y1 - y3) / (y1 - y2));
    const double phi = asin(sqrt((y1 - y2) / ((1 + z) * arg1 + y1)));
    const double phi0 = asin(sqrt((y1 - y2) / (arg1 + y1)));
    return (1. + z) / m_sqok * sin(g * (ellint_1(k, phi0) - ellint_1(k,
        phi)));
   }

   double flrw_nat_A2_1::age(double z) const
   {
     return 0;
   }

   double flrw_nat_A2_2::age(double z) const
   {
     return 0;
   }

  } // namespace impl

} // namespace milia

