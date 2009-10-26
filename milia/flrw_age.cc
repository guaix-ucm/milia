/*
 * Copyright 2008-2009 Sergio Pascual
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

// $Id$

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "flrw.h"
#include "flrw_prec.h"
#include "exception.h"

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>

using std::abs;

namespace
{
  double helper_fun_time(double z, void* pars)
  {
    milia::metrics::flrw* pmetric = static_cast<milia::metrics::flrw*> (pars);
    const double om = pmetric->get_matter();
    const double ol = pmetric->get_vacuum();
    return 1 / ((1 + z) * (sqrt(gsl_pow_2(1 + z) * (1 + om * z) - z * ol * (2
        + z))));
  }
}

namespace milia
{
  namespace metrics
  {

    double flrw::age() const
    {
      return m_uage;
    }

    double flrw::age(double z) const
    {
      switch (m_case)
      {
        case OM_OV_0:
          return m_t_h / (1 + z);
        case OV_1:
        case OV_2:
        case OV_EDS:
          return tolz(z);
        case OM:
          return tomz(z);
        case A1:
          return ta1(z);
        case A2_1:
        case A2_2:
          return ta2(z);
        case OM_OV_1:
          return tb(z);
        case OM_DS:
        case NO_CASE:
          break;
      }
      return -1;
    }

    // ov == 0 CASE: OV_1, OV_2, OV_EDS
    double flrw::tolz(double z) const
    {
      const double prez = sqrt(1.0 + m_om * z);
      switch (m_case)
      {
        case OV_1:
          //ov=0 0 < om < 1
          return m_t_h * (prez / (1.0 + z) - m_om / m_sqok * atanh(m_sqok
              / prez)) / m_ok;
        case OV_2:
          //ov=0 om > 1
          return m_t_h * (prez / (1.0 + z) - m_om / m_sqok
              * atan(m_sqok / prez)) / m_ok;
        case OV_EDS:
          // ov = 0 om = 1
          return m_t_h * 2.0 / (3.0 * (1 + z) * sqrt(1 + z));
        default:
          break;
      }
      return -1.;
    }

    // om=0 CASE: OM
    double flrw::tomz(double z) const
    {
      return m_t_h * asinh(1.0 / ((1 + z) * sqrt(1.0 / m_ov - 1.0))) / sqrt(
          m_ov);
    }

    // CASE A1
    double flrw::ta1(double z) const
    {
      gsl_error_handler_t* oldhandler = gsl_set_error_handler_off();
      int status = 0;
      gsl_sf_result result;
      //
      const double vk =
          cbrt(m_kap * (m_crit - 1) + sqrt(m_crit * (m_crit - 2)));
      const double y1 = (m_kap * (vk + 1 / vk) - 1) / 3.;
      const double A = sqrt(y1 * (3 * y1 + 2));
      // Parameters of the elliptical functions
      const double k = sqrt((2 * A + m_kap * (1 + 3 * y1)) / (4 * A));
      double arg0 = m_kap * y1 + m_om * (1 + z) / abs(m_ok);
      double phi = acos((arg0 - A) / (arg0 + A));
      // Selecting between cases
      // these conditions must hold
      // abs(k) <= 1 and n * gsl_pow_2(sin(phi_z)) < 1
      const double sin_phi = sin(phi);
      // Parameters of the elliptical function of third kind
      const double n_10 = gsl_pow_2(A + m_kap * y1) / (4 * A * m_kap * y1);
      const double n_8 = y1 * (1 + y1) / gsl_pow_2(A - m_kap * y1);

      double arg1 = (1 + z) * m_om / m_ok;

      const double crit8 = 1 - n_8 * gsl_pow_2(sin_phi);
      const double crit10 = 1 - n_10 * gsl_pow_2(sin_phi);

      if (abs(crit10) < FLRW_EQ_TOL)
      {
        // check if there's a node in eq 8 also, try eq 22 if so
        if (abs(crit8) < FLRW_EQ_TOL)
        {
          //  Equation 22, a very special case of b = 27*(2+sqrt(2))/8.
          const double phi = acos(-1 - arg1 / M_SQRT2 + 1 - arg1);
          const double arg2 = (1 - arg1) * sqrt(gsl_pow_2(arg1) + (1
              + M_SQRT1_2) * (arg1 + M_SQRT1_2));
          const double arg3 = (M_SQRT2 - 1) * (arg1 + M_SQRT2 + 1) * sqrt((1
              + M_SQRT1_2) * (M_SQRT1_2 - arg1));
          gsl_set_error_handler(oldhandler);
          return m_t_h / (4 * sqrt(m_ov)) * ((M_SQRT2 - 1) * gsl_sf_ellint_F(
              phi, 0.5 * sqrt(1 + 2 * M_SQRT2), ELLIP_PREC) + log(abs((arg2
              + arg3) / (arg2 - arg3))));
        }
        // if the critical parameter is negative, we have imaginary terms
        // for the moment we must integrate
        else if (crit8 < 0)
        {
          return ti(z);
        }
        // if not, go ahead with equation 8
        else
        {
          double arg2 = 2 * m_kap * y1 * sqrt((1 + y1) * ((1 + y1) * gsl_pow_2(
              y1) - gsl_pow_2(arg1) * (1 + arg1)));
          arg1 *= arg1 * (A - m_kap * y1) - 2 * m_kap * (1 + y1)
              * gsl_pow_2(y1);
          const double hm = arg1 - arg2;
          const double hp = arg1 + arg2;
          status = gsl_sf_ellint_P_e(phi, k, -n_8, ELLIP_PREC, &result);
          gsl_set_error_handler(oldhandler);
          if (status)
          {
            return ti(z);
          }
          arg2 = (A - m_kap) / (y1 * (1.0 + y1) * sqrt(A)) * result.val;
          const double arg3 = log(abs(hm / hp)) / (m_kap * y1 * sqrt(m_kap
              * (y1 + 1)));
          const double pre = 0.5 * m_om / sqrt(gsl_pow_3(abs(m_ok)));
          return m_t_h * pre * (arg1 + arg2 + arg3);
        }
      }
      else if (crit10 < 0)
      {
        return ti(z);
      }
      // if not, go ahead with equation 10
      else
      {
        // Equation 10
        const double hm = sqrt(((1 + y1) * (y1 - arg1)) / (gsl_pow_2(y1) + (1
            + arg1) * (y1 + arg1)));
        arg1 = -gsl_sf_ellint_F(phi, k, ELLIP_PREC) / (A + m_kap * y1);
        status = gsl_sf_ellint_P_e(phi, k, -n_10, ELLIP_PREC, &result);
        gsl_set_error_handler(oldhandler);
        if (status)
        {
          return ti(z);
        }
        const double arg2 = -0.5 * (A - m_kap * y1) / (m_kap * y1 * (A + m_kap
            * y1)) * result.val;
        const double arg3 = -0.5
            * (sqrt(A / (m_kap * (y1 + 1))) / (m_kap * y1)) * log(abs(
            (1.0 - hm) / (1.0 + hm)));
        return m_t_h * m_om / (abs(m_ok) * sqrt(A * abs(m_ok))) * (arg1 + arg2
            + arg3);
      }
      return -1.0;
    }

    double flrw::ta2(double z) const
    {
      // Equation 19, a very special case of b = 2.
      if (m_case == A2_1)
      {
        const double arg = (1 + z) * m_om / m_ok;
        return -m_t_h * (M_SQRT3 * log(1 - 2 / (sqrt(1. / 3 - arg) + 1)) + log(
            1 + 2. / (sqrt(1 - 3 * arg) - 1))) / sqrt(m_ov);
      }
      // Equation 15.
      if (m_case == A2_2)
      {
        int status = 0;
        gsl_sf_result result;
        const double arg = acos(1 - m_crit) / 3;
        const double arg1 = cos(arg) / 3;
        const double arg2 = sin(arg) / M_SQRT3;
        const double y1 = -1.0 / 3.0 + arg1 + arg2;
        const double y2 = -1.0 / 3.0 - 2 * arg1;
        const double arg3 = y1 - y2;
        const double y3 = y1 - 2 * arg2;
        const double phi = asin(sqrt(arg3 / (y1 - m_om * (1 + z) / m_ok)));
        const double k = sqrt((y1 - y3) / arg3);
        const double n = -y1 / arg3;
        const double arg4 = y1 * abs(m_ok) * sqrt(-arg3 * m_ok);
        status = gsl_sf_ellint_P_e(phi, k, n, ELLIP_PREC, &result);
        if (status)
        {
          return ti(z);
        }
        return m_t_h * 2.0 * m_om / arg4 * (result.val - gsl_sf_ellint_F(phi,
            k, ELLIP_PREC));
      }
      return -1.0;
    }

    double flrw::tb(double z) const
    {
      // om + ol == 1
      return m_t_h * 2 * asinh(sqrt(m_ov / (m_om * (1 + z * (3 + z * (3
          + z)))))) / (3 * sqrt(m_ov));
    }

    double flrw::ti(double z) const
    {
      gsl_error_handler_t* oldhandler = gsl_set_error_handler_off();
      gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
      double result, error;
      gsl_function F;
      F.function = &helper_fun_time;
      F.params = static_cast<void*> (const_cast<flrw*> (this));
      const int status = gsl_integration_qagiu(&F, z, 0, 1e-7, 1000, w,
          &result, &error);
      gsl_integration_workspace_free(w);
      gsl_set_error_handler(oldhandler);
      if (status)
      {
        throw milia::exception(std::string("gsl error: ")
            + gsl_strerror(status));
      }
      return m_t_h * result;
    }

  } // namespace metrics
} // namespace milia
