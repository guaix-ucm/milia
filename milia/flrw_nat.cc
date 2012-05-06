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

#include "flrw_nat.h"
#include "exception.h"
#include "flrw_prec.h"

#include <cmath>
#include <sstream>
#include <limits>

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

    flrw_nat::flrw_nat(double m, double v) :
          m_om(m), m_ov(v), m_ok(1 - m_om - m_ov), m_sqok(sqrt(abs(m_ok))),
          m_kap(m_ok > 0 ? -1 : 1), m_crit(-13.5 * gsl_pow_2(m_om) * m_ov
              / gsl_pow_3(m_ok)), m_case(select_case(m_om, m_ov, m_ok, m_crit))
    {
      gsl_set_error_handler_off();

      // om < 0 not allowed
      if (m_om < -FLRW_EQ_TOL)
        throw milia::exception("Matter density < 0 not allowed");

      //ov < 0 makes the universe recollapse
      if (m_ov < -FLRW_EQ_TOL)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      if (m_ov >= 1 && m_crit <= 2)
        throw milia::no_big_bang("No Big Bang"); // No Big bang with these parameters

      m_flags.time_ends = false;
      m_flags.time_begins = true;

      switch(m_case) {
      case OM_DS:
        m_flags.time_begins = false;
        m_flags.time_begin_scale = std::numeric_limits<double>::infinity();
      default:
        m_flags.time_begin_scale = age(0);
      }

    }

    std::string flrw_nat::to_string() const
    {
      std::stringstream out;
      out << "flrw_nat(matter=" << m_om << ", vacuum=" << m_ov << ")";
      return out.str();
    }

    bool flrw_nat::does_recollapse(double matter, double vacuum)
    {
      if (vacuum < 0)
        return true;
      if (matter < 1)
        return false;
      const double critical = 4 * matter * gsl_pow_3(cos(1. / 3. * acos(1.
          / matter - 1.) + 4 * M_PI / 3.));
      if (vacuum > critical)
        return false;
      return true;
    }

    void flrw_nat::set_matter(double M)
    {
      if (M < -FLRW_EQ_TOL)
      {
        throw milia::exception("Matter density < 0 not allowed");
      }

      const double OK = 1 - m_ov - M;
      double B = -13.5 * gsl_pow_2(M) * m_ov / gsl_pow_3(OK);

      if (M >= 1 && B <= 2)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      if (m_ov >= 1 && B <= 2)
        throw milia::no_big_bang("No Big Bang"); // No Big bang with this parameters

      m_om = M;
      m_ok = OK;
      m_crit = B;
      m_kap = (m_ok > 0 ? -1 : 1);
      m_sqok = sqrt(abs(m_ok));
      m_case = select_case(m_om, m_ov, m_ok, m_crit);
      m_flags.time_begin_scale = age(0);
    }

    void flrw_nat::set_vacuum(double L)
    {
      if (L < -FLRW_EQ_TOL)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      const double OK = 1 - m_om - L;
      const double B = -13.5 * gsl_pow_2(m_om) * L / gsl_pow_3(OK);

      if (m_om >= 1 && B <= 2)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      if (L >= 1 && B <= 2)
        throw milia::no_big_bang("No Big Bang"); // No Big bang with this parameters

      m_ov = L;
      m_ok = OK;
      m_crit = B;
      m_kap = (m_ok > 0 ? -1 : 1);
      m_kap = (m_ok > 0 ? -1 : 1);
      m_sqok = sqrt(abs(m_ok));
      m_case = select_case(m_om, m_ov, m_ok, m_crit);
      m_flags.time_begin_scale = age(0);
    }

    flrw_nat::ComputationCases flrw_nat::select_case(double om, double ov, double ok,
        double crit)
    {
      const bool l3 = (abs(om) < FLRW_EQ_TOL);
      const bool l4 = (abs(ov) < FLRW_EQ_TOL);
      // om = ov = 0
      if (l3 && l4)
        return OM_OV_0;
      // ov = 0 om == 1. Einstein-de Sitter Universe
      if (l4 && abs(om - 1) < FLRW_EQ_TOL)
        return OV_EDS;
      // ov=0 0<om<1
      if (l4 && om < 1)
        return OV_1;
      // ov=0 om>1
      if (l4 && om > 1)
        return OV_2;
      // om = 0 ov = 1 de Sitter Universe
      if (l3 and (abs(ov - 1) < FLRW_EQ_TOL))
        return OM_DS;
      // om = 0 0 < ov < 1
      if (l3 and (ov > 0) and (ov < 1))
        return OM;
      // om + ov == 1, flat Universe
      if (abs(ok) < FLRW_EQ_TOL)
        return OM_OV_1;
      // crit == 2
      if (abs(crit - 2) < FLRW_EQ_TOL)
        return A2_1;
      // om+ov != 1 and crit < 0 or crit > 2
      if ((crit < 0) || (crit > 2))
        return A1;

      // om+ov != 1 and 0 < crit < 2
      if ((crit > 0) && (crit < 2))
        return A2_2;
      //
      return NO_CASE;
    }

    double flrw_nat::hubble(double z) const
    {
      return sqrt(m_ov * gsl_pow_3(1 + z) + m_ok * gsl_pow_2(1 + z)
          + m_ov);
    }

    double flrw_nat::lt(double z) const
    {
      // de Sitter's Universe age is infinity
      // but look-back time is valid
      switch (m_case) {
       case OM_DS:
        return log(1 + z);
       default:
        return m_flags.time_begin_scale - age(z);
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
        {
          const double r = dm;
          return asinc(m_kap, m_sqok, r);
        }
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
          return gsl_pow_3(dm) / 3.0;
        default:
        {
          const double r = dm;
          return (r * sqrt(1 + m_ok * gsl_pow_2(r)) - asinc(
              m_kap, m_sqok, r)) / (2 * m_ok);
        }

      }
    }

    double flrw_nat::sinc(double k, double a, double x)
    {
      if (k > 0)
        return sin(a * x) / a;
      if (k < 0)
        return sinh(a * x) / a;
      return -1;
    }

    double flrw_nat::asinc(double k, double a, double x)
    {
      if (k > 0)
        return asin(a * x) / a;
      if (k < 0)
        return asinh(a * x) / a;
      return -1;
    }

  } //namespace metrics
} //namespace milia

std::ostream& operator<<(std::ostream& os, milia::metrics::flrw_nat& iflrw)
{
  os << iflrw.to_string();
  return os;
}
