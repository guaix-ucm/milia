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

#include "flrw_nat_impl.h"
#include "util.h"

using boost::math::pow;

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

  } // namespace impl

} // namespace milia

