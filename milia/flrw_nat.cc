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
#include "metric.h"
#include "util.h"

#include "flatmodel.h"

using std::abs;
using boost::math::asinh;
using boost::math::pow;

namespace milia
{

    flrw_nat::flrw_nat(double m, double v) :
      m_om(m), m_ov(v), m_ok(1 - m_om - m_ov), m_sqok(sqrt(abs(m_ok)))
    {

      //om < 0 not allowed
      if (m_om < 0)
        throw std::domain_error("Matter density < 0 not allowed");

      //ov < 0 not allowed
      if (m_ov < 0)
        throw std::domain_error("Vacuum density < 0 not allowed");

      m_crit = -13.5 * pow<2> (m_om) * m_ov / (pow<3> (m_ok));

      m_kap = m_ok > 0 ? -1 : 1;
      m_case = select_case();

      /*
          NO_CASE, // error condition
          OM_OV_0, //om = ov = 0
          OV_1, //ov = 0 0 < om < 1
          OV_2, //ov = 0 om > 1
          OV_EDS, //ov = 0 om = 1, Einstein-de Sitter Universe
          OM, //om = 0 0 < ov < 1
          OM_DS, //om = 0 ov = 1, de Sitter Universe
          OM_OV_1, //om + ol = 1
          A1, //om+ov != 1 b < 0 || b > 2
          A2_1, //om+ov != 1 b = 2 Lower b=2 curve
          A2_2, //om+ov != 1 0 < b < 2 Only valid in the lower part
      */

      if (m_case == A2_1 or m_case == A2_2)
      {
        if (m_ov > m_om) // we are in the upper region
          throw std::domain_error("No Big Bang"); // No Big bang with this parameters

      }

      m_uage = m_case != OM_DS ? age(0) : 0;
    }

    std::string flrw_nat::to_string() const
    {
      std::stringstream out;
      out << "flrw_nat(matter=" << m_om << ", vacuum="
          << m_ov << ")";
      return out.str();
    }



    bool flrw_nat::does_recollapse(double matter, double vacuum)
    {
       return check_recollapse(matter, vacuum);
    }

    void flrw_nat::set_matter(double matter)
    {
      if (matter < 0)
      {
        throw std::domain_error("Matter density < 0 not allowed");
      }

      const double OK = 1 - m_ov - matter;
      double B = -13.5 * pow<2> (matter) * m_ov / pow<3> (OK);

      if (matter >= 1 && B <= 2)
        throw std::domain_error("The Universe recollapses"); // Recollapse

      if (m_ov >= 1 && B <= 2)
        throw std::domain_error("No Big Bang"); // No Big bang with this parameters

      // Time and space scale
      m_om = matter;
      m_ok = OK;
      m_sqok = sqrt(abs(OK));
      m_crit = B;
      m_kap = (m_ok > 0 ? -1 : 1);
      m_case = select_case();
      m_uage = m_case != OM_DS ? age() : 0;
    }

    void flrw_nat::set_vacuum(double vacuum)
    {
      if (vacuum < 0)
        throw std::domain_error("The Universe recollapses"); // Recollapse

      const double OK = 1 - m_om - vacuum;
      double B = -13.5 * pow<2> (m_om) * vacuum / pow<3> (OK);

      if (m_om >= 1 && B <= 2)
        throw std::domain_error("The Universe recollapses"); // Recollapse

      if (vacuum >= 1 && B <= 2)
        throw std::domain_error("No Big Bang"); // No Big bang with this parameters

      m_ov = vacuum;
      m_ok = OK;
      m_sqok = sqrt(abs(OK));
      m_crit = B;
      m_kap = (m_ok > 0 ? -1 : 1);
      m_case = select_case();
      m_uage = m_case != OM_DS ? age(0) : 0;
    }

    flrw_nat::ComputationCases flrw_nat::select_case() const
    {
      const bool l3 = (abs(m_om) < FLRW_EQ_TOL);
      const bool l4 = (abs(m_ov) < FLRW_EQ_TOL);
      // om = ov = 0 Degenerate case
      if (l3 && l4)
        return OM_OV_0;
      // ov = 0 om == 1. Einstein-de Sitter Universe
      if (l4 && abs(m_om - 1) < FLRW_EQ_TOL)
        return OV_EDS;
      // ov=0 0<om<1
      if (l4 && m_om < 1)
        return OV_1;
      // ov=0 om>1
      if (l4 && m_om > 1)
        return OV_2;

      // om = 0 ov = 1 de Sitter Universe
      if (l3 and (abs(m_ov - 1) < FLRW_EQ_TOL))
        return OM_DS;
      // om = 0 0 < ov < 1
      if (l3 and (m_ov > 0) and (m_ov < 1))
        return OM;
      // om + ov == 1, flat Universe
      if (abs(m_ok) < FLRW_EQ_TOL)
        return OM_OV_1;
      // b == 2
      if (abs(m_crit - 2) < FLRW_EQ_TOL)
        return A2_1;
      // om+ov != 1 b<0 || b>2
      if ((m_crit < 0) || (m_crit > 2))
        return A1;

      // om+ov != 1 0 < b < 2
      if ((m_crit > 0) && (m_crit < 2))
        return A2_2;
      //
      return NO_CASE;
    }

    double flrw_nat::get_hubble(double z) const
    {
      return sqrt(m_om * pow<3> (1 + z) + m_ok * pow<2> (1 + z) + m_ov);
    }

    double flrw_nat::lt(double z) const
    {
      // de Sitter's Universe age is infinity
      // but look-back time is valid
      switch(m_case) {
        case OM_DS:
         return log(1 + z);
        default:
         return m_uage - age(z);
      }
    }

    double flrw_nat::dc(double z, double dm) const
    {
      switch (m_case)
      {
        // Flat cases
        case OV_EDS: // EdS
        case OM_DS: // dS
        case OM_OV_1: // Flat
          return dm;
        default:
          return asinc(m_kap, m_sqok, dm);
      }
    }

    double flrw_nat::vol(double z, double dm) const
    {
      switch (m_case)
      {
        // Flat cases
        case OV_EDS: // EdS
        case OM_DS: // dS
        case OM_OV_1: // Flat
          return pow<3> (dm) / 3.0;
        default:
          return (dm * sqrt(1 + m_ok * pow<2> (dm)) - asinc(m_kap,
              m_sqok, dm)) / (2 * m_ok);
      }
    }

} //namespace milia

std::ostream& operator<<(std::ostream& os, milia::flrw_nat& iflrw)
{
  os << iflrw.to_string();
  return os;
}

namespace milia {
    flrw_nat_new::flrw_nat_new(double m, double v) {
    
      m_impl.reset(impl::flrw_nat_impl::construct(m, v));
    }

}
