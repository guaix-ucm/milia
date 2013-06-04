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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>
#include <iostream>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp>

#include "flatmodel.h"

using boost::math::asinh;
using boost::math::pow;
using boost::math::ellint_1;
using boost::math::cbrt;

namespace
{
  const double M_SQRT3 = sqrt(3);
  const double M_4THRT3 = sqrt(M_SQRT3);
}

namespace milia
{

  namespace impl {
    const char* flrw_nat_OV_EDS::model() const {
      return "OV_EDS";
    }

    double flrw_nat_OV_EDS::dl(double z) const
    {
      return 2 * (1 + z - std::sqrt(1 + z));
    }

    double flrw_nat_OV_EDS::age(double z) const
    {
      return 2 / (3 * (1 + z) * std::sqrt(1 + z));
    }

    const char* flrw_nat_OM_DS::model() const {
      return "OM_DS";
    }

    double flrw_nat_OM_DS::dl(double z) const
    {
      return z * (z + 1);
    }

    double flrw_nat_OM_DS::age(double z) const
    {
      return 0;
    }

    double flrw_nat_OM_DS::lt(double z) const
    {
      return std::log(1 + z);
    }

    const char* flrw_nat_OM_OV_1::model() const {
      return "OM_OV_1";
    }

    double flrw_nat_OM_OV_1::age(double z) const
    {
      return 2. / (3. * std::sqrt(m_ov)) * asinh(std::sqrt((1 / m_om - 1) / pow<3>(1 + z)));
    }

    double flrw_nat_OM_OV_1::dl(double z) const
    {
      const double k = sqrt(0.5 + 0.25 * M_SQRT3);
      const double c1 = M_4THRT3;
      const double arg0 = cbrt((1 / m_om - 1));
      const double down = 1 + (1 + M_SQRT3) * arg0;
      const double up = 1 + (1 - M_SQRT3) * arg0;
      const double phi = acos((z + up) / (z + down));
      const double phi0 = acos(up / down);
      return (1 + z) / (c1 * sqrt(m_om) * sqrt(arg0)) * (ellint_1(k, phi0)
           - ellint_1(k, phi));
    }
  } //namespace impl

} //namespace milia
