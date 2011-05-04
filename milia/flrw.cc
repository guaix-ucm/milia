/*
 * Copyright 2008-2011 Sergio Pascual
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

#include "flrw.h"
#include "exception.h"
#include "flrw_prec.h"

#include <cmath>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#ifndef HAVE_ASINH
#define asinh gsl_asinh
#endif

#ifndef HAVE_ATANH
#define atanh gsl_atanh
#endif

namespace milia
{
  namespace metrics
  {
    using std::abs;

    flrw::flrw(double h, double m, double v) :
      flrw_nat(m, v),
      m_hu(h), m_r_h(ms_hubble_radius / m_hu), m_t_h(ms_hubble_time / m_hu)
    {
      if (m_hu <= 0)
        throw milia::exception("Hubble constant <= 0 not allowed");
    }

    std::string flrw::to_string() const
    {
      std::stringstream out;
      out << "flrw(hubble=" << m_hu << ", matter=" << get_matter()
          << ", vacuum=" << get_vacuum() << ")";
      return out.str();
    }

    void flrw::set_hubble(double hubble)
    {
      if (hubble > 0)
      {
        m_hu = hubble;
        m_r_h = ms_hubble_radius / m_hu;
        m_t_h = ms_hubble_time / m_hu;
        return;
      }
      else
      {
        throw milia::exception("Hubble constant <= 0 not allowed");
      }
    }

    double flrw::DM(double z) const
    {
      return 5 * log10(dl(z)) + 25;
    }

    double flrw::angular_scale(double z) const
    {
      // 206264.8062 converts arcsconds to radians
      const double arcsec_to_rad = 206264.8062;
      return da(z) * 1e6 / arcsec_to_rad;
    }

  } //namespace metrics
} //namespace milia

std::ostream& operator<<(std::ostream& os, milia::metrics::flrw& iflrw)
{
  os << iflrw.to_string();
  return os;
}
