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

#ifndef MILIA_FLRW_NAT_FLATMODEL_H
#define MILIA_FLRW_NAT_FLATMODEL_H

#include "flrw_nat_impl.h"

namespace milia
{
  namespace impl {

    class flrw_nat_flat: public milia::impl::flrw_nat_impl
    {
     public:
       flrw_nat_flat(double matter, double vacuum) :
         flrw_nat_impl(matter, vacuum)
       {}

        double dc(double z) const {
          return this->dm(z);
        }

        double vol(double z) const {
          const double lm = this->dm(z);
          return lm * lm * lm / 3.0;
        }
    };

    class flrw_nat_OV_EDS : public flrw_nat_flat
    {
      public:
        flrw_nat_OV_EDS() :
         flrw_nat_flat(1.0, 0.0)
        {}


        double dl(double z) const;
        double age(double z) const;
    };
 
  } // namespace impl

} // namespace milia

#endif /* MILIA_FLRW_NAT_FLATMODEL_H */
