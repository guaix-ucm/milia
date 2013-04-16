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

#include "flatmodel.h"

namespace milia
{

  namespace impl {

    double flrw_nat_OV_EDS::dl(double z) const
    {
      return 2 * (1 + z - std::sqrt(1 + z));
    }

    double flrw_nat_OV_EDS::age(double z) const
    {
      return 2 / (3 * (1 + z) * std::sqrt(1 + z));
    }

  } //namespace impl

} //namespace milia
