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

#include <cmath>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>

#include "flrw_nat.h"
#include "flrw_prec.h"
#include "exception.h"

using std::sqrt;
using std::atan;
using std::log;
using std::abs;

using boost::math::ellint_1;
using boost::math::ellint_3;
using boost::math::asinh;
using boost::math::atanh;
using boost::math::cbrt;
using boost::math::pow;

namespace
{
  //const double M_SQRT3 = sqrt(3);
  const double M_4THRT3 = sqrt(M_SQRT3);
  double helper_fun_time(double z, void* pars)
  {
    milia::metrics::flrw_nat* pmetric =
        static_cast<milia::metrics::flrw_nat*> (pars);
    const double om = pmetric->get_matter();
    const double ol = pmetric->get_vacuum();
    return 1. / ((1 + z) * (sqrt(pow<2> (1 + z) * (1 + om * z) - z * ol * (2
        + z))));
  }
}

namespace milia
{
  namespace metrics
  {
    double flrw_nat::age() const
    {
      return m_uage;
    }

    double flrw_nat::age(double z) const
    {
      switch (m_case)
      {
        case OM_OV_0:
          return 1.0 / (1 + z);
          break;
        case OV_1:
        case OV_2:
        case OV_EDS:
          return tolz(z);
          break;
        case OM:
          return tomz(z);
          break;
        case A1:
          return ta1(z);
          break;
        case A2_1:
        case A2_2:
          return ta2(z);
          break;
        case OM_OV_1:
          return tb(z);
          break;
      }
      return -1;
    }

    // ol = 0 CASE: OV_1, OV_2, OV_EDS
    double flrw_nat::tolz(double z) const
    {
      const double pre0 = 1 - m_om;
      const double prez = sqrt(1 + m_om * z);
      switch (m_case)
      {
        case OV_1:
          //ov = 0 0 < om < 1
          return (prez / (1 + z) - m_om / sqrt(pre0) * atanh(sqrt(pre0) / prez))
              / pre0;
        case OV_2:
          //ol = 0 om > 1
          return (prez / (1 + z) - m_om / sqrt(-pre0)
              * atan(sqrt(-pre0) / prez)) / pre0;
        case OV_EDS:
          // ol = 0 om = 1
          return 2 / (3 * (1 + z) * sqrt(1 + z));
      }
      return -1.;
    }

    // om=0 CASE: OM
    double flrw_nat::tomz(double z) const
    {
      return asinh(1 / ((1 + z) * sqrt(1 / m_ov - 1))) / sqrt(m_ov);
    }

    // CASE A1
    double flrw_nat::ta1(double z) const
    {
      const double vk = cbrt(m_kap * (m_crit - 1) + sqrt(m_crit * (m_crit - 2)));
      const double y1 = (m_kap * (vk + 1 / vk) - 1) / 3.;
      const double A = sqrt(y1 * (3 * y1 + 2));
      // Parameters of the elliptical functions
      const double k = sqrt((2 * A + m_kap * (1 + 3 * y1)) / (4 * A));
      double arg0 = m_kap * y1 + m_om * (1 + z) / abs(m_ok);
      double phi = acos((arg0 - A) / (arg0 + A));

      // Selecting between cases
      // these conditions must hold
      // abs(k) <= 1 and n * pow<2>(sin(phi_z)) < 1

      const double sin_phi = sin(phi);
      const double n_10 = pow<2> (A + m_kap * y1) / (4 * A * m_kap * y1);
      const double n_8 = y1 * (1 + y1) / pow<2> (A - m_kap * y1);

      double arg2, arg3;
      double hm, hp;
      double pre;
      double arg1 = (1 + z) * m_om / m_ok;

      const double crit8 = 1 - n_8 * pow<2> (sin_phi);
      const double crit10 = 1 - n_10 * pow<2> (sin_phi);
      // If there's a node in eq 10, try eq 8
      if (abs(crit10) < FLRW_EQ_TOL)
      {
        // check if there's a node in eq 8 also, try eq 22 if so
        if (abs(crit8) < FLRW_EQ_TOL)
        {
          //  Equation 22, a very special case of b = 27*(2+sqrt(2))/8.
          phi = acos(-1 - arg1 / M_SQRT2 + 1 - arg1);
          //      std::cout << "eq 22" << std::endl;
          const double arg2 = (1 - arg1) * sqrt(arg1 * arg1
              + (arg1 + M_SQRT1_2) * (1 + M_SQRT1_2));
          const double arg3 = (M_SQRT2 - 1) * (arg1 + M_SQRT2 + 1) * sqrt((1
              + M_SQRT1_2) * (M_SQRT1_2 - arg1));
          return 0.25 / sqrt(m_ov) * ((M_SQRT2 - 1) * ellint_1(0.5 * sqrt(1 + 2
              * M_SQRT2), phi) + log(abs((arg2 + arg3) / (arg2 - arg3))));
        }
        //        EQUATION 8.
        // if the critical parameter is negative, we have imaginary terms
        // for the moment we must integrate
        else if (crit8 < 0)
        {
          return ti(z);

        }
        // if not, go ahead with equation 8
        else
        {
          //       std::cout << "eq 8" << std::endl;
          arg2 = arg1 * arg1 * (1 + arg1);
          arg2 = sqrt((1 + y1) * ((1 + y1) * pow<2> (y1) - arg2));
          arg2 *= 2 * m_kap * y1;
          arg1 *= arg1 * (A - m_kap * y1) - 2 * m_kap * (1 + y1) * pow<2> (y1);
          hm = arg1 - arg2;
          hp = arg1 + arg2;
          arg1 = ellint_1(k, phi) / (m_kap * y1 * sqrt(A));
          arg2 = (A - m_kap) / (y1 * (1 + y1) * sqrt(A))
              * ellint_3(k, n_8, phi);
          arg3 = log(abs(hm / hp)) / (m_kap * y1 * sqrt(m_kap * (y1 + 1)));
          pre = 0.5 * m_om / (abs(m_ok) * m_sqok);
          return pre * (arg1 + arg2 + arg3);
        }
      }
      else if (crit10 < 0)
      {
        return ti(z);
      }

      else
      {
        // Equation 10
        //        std::cout << "eq 10" << std::endl;
        hm = sqrt(((1 + y1) * (y1 - arg1)) / (pow<2> (y1) + (1 + arg1) * (y1
            + arg1)));
        arg1 = -ellint_1(k, phi) / (A + m_kap * y1);
        arg2 = -0.5 * (A - m_kap * y1) / (m_kap * y1 * (A + m_kap * y1))
            * ellint_3(k, n_10, phi);
        arg3 = -0.5 * (sqrt(A / (m_kap * (y1 + 1))) / (m_kap * y1)) * log(abs(
            (1.0 - hm) / (1.0 + hm)));
        return m_om / (sqrt(A * abs(pow<3> (m_ok)))) * (arg1 + arg2 + arg3);
      }
      return -1.0;
    }

    double flrw_nat::ta2(double z) const
    {
      //       EQUATION 19, a very special case of b = 2.
      if (m_case == A2_1)
      {
        const double arg = (1 + z) * m_om / m_ok;
        return -(M_SQRT3 * log(1 - 2 / (sqrt(1. / 3. - arg) + 1)) + log(1 + 2
            / (sqrt(1 - 3 * arg) - 1))) / sqrt(m_ov);
      }
      //       EQUATION 15.
      if (m_case == A2_2)
      {
        const double arg = acos(1 - m_crit) / 3;
        const double arg1 = cos(arg) / 3;
        const double arg2 = sin(arg) / M_SQRT3;
        const double y1 = -1. / 3. + arg1 + arg2;
        const double y2 = -1. / 3. - 2 * arg1;
        const double arg3 = y1 - y2;
        const double y3 = y1 - 2 * arg2;
        const double phi = asin(sqrt(arg3 / (y1 - m_om * (1 + z) / m_ok)));
        const double k = sqrt((y1 - y3) / arg3);
        const double n = -y1 / arg3;
        const double arg4 = y1 * abs(m_ok) * sqrt(-arg3 * m_ok);
        return 2 * m_om / arg4 * (ellint_3(k, n, phi) - ellint_1(k, phi));
      }
      return -1.0;
    }

    double flrw_nat::tb(double z) const
    {
      // om + ov = 1
      return (2. / (3. * sqrt(1 - m_om))) * asinh(sqrt((1. / m_om - 1)
          / (pow<3> (1 + z))));
    }
    double flrw_nat::ti(double z) const
    {
      gsl_error_handler_t* oldhandler = gsl_set_error_handler_off();
      gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
      double result, error;
      gsl_function F;
      F.function = &helper_fun_time;
      F.params = static_cast<void*> (const_cast<flrw_nat*> (this));
      const int status = gsl_integration_qagiu(&F, z, 0, 1e-7, 1000, w,
          &result, &error);
      gsl_integration_workspace_free(w);
      gsl_set_error_handler(oldhandler);
      if (status)
      {
        throw milia::exception(std::string("gsl error: ")
            + gsl_strerror(status));
      }
      return result;
    }
  } // namespace metrics
} // namespace milia
