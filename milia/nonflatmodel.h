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

        double dc(double z) const;

        double vol(double z) const;
    };

    class flrw_nat_OM_OV_0: public flrw_nat_nonflat 
    {
     public:
        flrw_nat_OM_OV_0(): flrw_nat_nonflat(0.0, 0.0)
        {}

        double dl(double z) const;

        double age(double z) const;

    };

  class flrw_nat_OV: public flrw_nat_nonflat 
  {
    public:
      flrw_nat_OV(double matter) : flrw_nat_nonflat(matter, 0.0)
      {}

      double dl(double z) const;
  };

  class flrw_nat_OV_1: public flrw_nat_OV 
  {
    public:
      flrw_nat_OV_1(double matter) : flrw_nat_OV(matter)
      {}

      double age(double z) const;
  };

  class flrw_nat_OV_2: public flrw_nat_OV 
  {
    public:
      flrw_nat_OV_2(double matter) : flrw_nat_OV(matter)
      {}

      double age(double z) const;
  };

  class flrw_nat_OM: public flrw_nat_nonflat 
  {
    public:
      flrw_nat_OM(double vacuum) : flrw_nat_nonflat(0.0, vacuum)
      {}

      double dl(double z) const;
      double age(double z) const;
  };

  class flrw_nat_A1: public flrw_nat_nonflat 
  {
    public:
      flrw_nat_A1(double matter, double vacuum) : flrw_nat_nonflat(matter, vacuum)
      {}

      double dl(double z) const;
      double age(double z) const;
    private:
      double ti(double z) const;
  };

  class flrw_nat_A2: public flrw_nat_nonflat 
  {
    public:
      flrw_nat_A2(double matter, double vacuum) : flrw_nat_nonflat(matter, vacuum)
      {}

      double dl(double z) const;
  };

  class flrw_nat_A2_1: public flrw_nat_A2
  {
    public:
      flrw_nat_A2_1(double matter, double vacuum) : flrw_nat_A2(matter, vacuum)
      {}

      double age(double z) const;
  };

  class flrw_nat_A2_2: public flrw_nat_A2
  {
    public:
      flrw_nat_A2_2(double matter, double vacuum) : flrw_nat_A2(matter, vacuum)
      {}

      double age(double z) const;
  };

  } // namespace impl

} // namespace milia

#endif /* MILIA_FLRW_NAT_NONFLATMODEL_H */
