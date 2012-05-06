/*
 * Copyright 2012 Sergio Pascual
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

#include <cmath>
#include <boost/math/special_functions/pow.hpp>

using boost::math::pow;

namespace milia
{
  bool check_recollapse(double matter, double vacuum) {

      if (vacuum < 0)
          return true;

      if ((vacuum > 0) and (matter > 1)) {

          const double crit = -13.5 * pow<2> (matter) * vacuum / pow<3> (1 - matter - vacuum);

          if ((crit > 0) and (crit <= 2)) {
             // critical region
             return true;
          }

      }
      return false;

  }

  bool check_bigbang(double matter, double vacuum) {
      if ((vacuum > 1) and (matter > 0)) {

          const double crit = -13.5 * pow<2> (matter) * vacuum / pow<3> (1 - matter - vacuum);

          if ((crit > 0) and (crit <= 2)) {
             // critical region
             return false;
          }

      }
      return true;
  }
  
} // namespace milia
