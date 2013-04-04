/*
 * Copyright 2008-2012 Sergio Pascual
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
#include <sstream>
#include <stdexcept>

#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/pow.hpp>

#include "flrw_nat.h"
#include "flrw_prec.h"

using std::abs;
using boost::math::asinh;
using boost::math::pow;

namespace milia
{
    double sinc(int k, double a, double x)
    {
      switch(k) {
      case 1:
        return sin(a * x) / a;
      case -1:
        return sinh(a * x) / a;
      case 0:
        return x;
      default:
        return 0;
      }
    }

    double asinc(int k, double a, double x)
    {
      switch(k) {
      case 1:
        return asin(a * x) / a;
      case -1:
        return asinh(a * x) / a;
      case 0:
        return x;
      default:
        return 0;
        }
    }

} //namespace milia
