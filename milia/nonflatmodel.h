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

#ifndef MILIA_FLRW_NAT_NONFLATMODEL_H
#define MILIA_FLRW_NAT_NONFLATMODEL_H

#include <boost/math/special_functions/pow.hpp>

#include "flrw_nat_impl.h"
#include "util.h"

using boost::math::pow;

namespace milia
{
  namespace impl {

    class flrw_nat_nonflat: public milia::impl::flrw_nat_impl
    {
     public:
       flrw_nat_nonflat(double matter, double vacuum) :
         flrw_nat_impl(matter, vacuum)
       {}

        double dc(double z) const {
          return asinc(m_kap, m_sqok, dm(z));
        }

        double vol(double z) const {
          const double lm = dm(z);
          return (lm * sqrt(1 + m_ok * pow<2> (lm)) - asinc(m_kap, m_sqok, lm)) / (2 * m_ok);
        }
    };

    class flrw_nat_OM_OV_0: public flrw_nat_nonflat 
    {
     public:
        flrw_nat_OM_OV_0(): flrw_nat_nonflat(0.0, 0.0)
        {}

        double dl(double z) const 
        {
          return 0.5 * z * (z + 2);
        }

        double age(double z) const
        {
          return 1.0 / (1 + z);
        }

    };


  } // namespace impl

} // namespace milia

#endif /* MILIA_FLRW_NAT_NONFLATMODEL_H */
