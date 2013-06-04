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

#include "nonflatmodel.h"
#include "util.h"
#include "flrw_prec.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/cbrt.hpp>

#include "flrw_nat_impl.h"
#include "util.h"

using std::sqrt;
using std::abs;
using std::sin;
using std::sinh;
using std::acos;
using std::atan;
using std::log;

using boost::math::ellint_1;
using boost::math::ellint_3;
using boost::math::cbrt;
using boost::math::asinh;
using boost::math::atanh;
using boost::math::pow;

namespace
{
  const double M_4THRT3 = std::sqrt(M_SQRT3);

  double helper_fun_time(double z, void* pars)
  {
    milia::impl::flrw_nat_impl* pmetric =
        static_cast<milia::impl::flrw_nat_impl*> (pars);
    const double om = pmetric->get_matter();
    const double ol = pmetric->get_vacuum();
    return 1. / ((1 + z) * (sqrt(pow<2> (1 + z) * (1 + om * z) - z * ol * (2
        + z))));
  }

}

namespace milia
{
  namespace impl {

   double flrw_nat_nonflat::dc(double z) const {
     return asinc(m_kap, m_sqok, dm(z));
   }

   double flrw_nat_nonflat::vol(double z) const {
     const double lm = dm(z);
     return (lm * sqrt(1 + m_ok * pow<2> (lm)) - asinc(m_kap, m_sqok, lm)) / (2 * m_ok);
   }

   double flrw_nat_OM_OV_0::dl(double z) const 
   {
     return 0.5 * z * (z + 2);
   }

   double flrw_nat_OM_OV_0::age(double z) const
   {
     return 1.0 / (1 + z);
   }

   double flrw_nat_OV::dl(double z) const
   {
     return 2 * ((2 - m_om * (1 - z) - (2 - m_om) * sqrt(1 + m_om * z)))
                   / pow<2> (m_om);
    }

   double flrw_nat_OV_1::age(double z) const 
   {
     const double pre0 = 1 - m_om;
     const double prez = sqrt(1 + m_om * z);
     return (prez / (1 + z) - m_om / sqrt(pre0) * atanh(sqrt(pre0) / prez)) / pre0;
   }

   double flrw_nat_OV_2::age(double z) const 
   {
     const double pre0 = 1 - m_om;
     const double prez = sqrt(1 + m_om * z);
     return (prez / (1 + z) - m_om / sqrt(-pre0) * atan(sqrt(-pre0) / prez)) / pre0;
   }
  
   double flrw_nat_OM::dl(double z) const 
   {
     return ((1 + z) / m_ov) * (1 + z - sqrt(m_ov + (1 - m_ov)
                 * pow<2> (1 + z)));
   }

   double flrw_nat_OM::age(double z) const 
   {
     return asinh(1 / ((1 + z) * sqrt(1 / m_ov - 1))) / sqrt(m_ov);
   }

   double flrw_nat_A1::dl(double z) const
   {
      const double v = cbrt(m_kap * (m_crit - 1) + sqrt(m_crit * (m_crit - 2)));
      const double y = (-1 + m_kap * (v + 1. / v)) / 3.;
      const double A = sqrt(y * (3 * y + 2));
      const double g = 1. / sqrt(A);
      const double k = sqrt(0.5 + 0.25 * pow<2> (g) * (v + 1. / v));
      const double sup = m_om / abs(m_ok);
      const double phi = acos(((1 + z) * sup + m_kap * y - A) / ((1 + z)
          * sup + m_kap * y + A));
      const double phi0 = acos((sup + m_kap * y - A)
          / (sup + m_kap * y + A));
      return (1 + z) / m_sqok * sinc(m_kap, 1.0, g * (ellint_1(k,
          phi0) - ellint_1(k, phi)));
   }

   double flrw_nat_A1::ti(double z) const
   {
      gsl_error_handler_t* oldhandler = gsl_set_error_handler_off();
      gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
      double result, error;
      gsl_function F;
      F.function = &helper_fun_time;
      F.params = static_cast<void*> (const_cast<flrw_nat_A1*> (this));
      const int status = gsl_integration_qagiu(&F, z, 0, 1e-7, 1000, w,
          &result, &error);
      gsl_integration_workspace_free(w);
      gsl_set_error_handler(oldhandler);
      if (status)
      {
        throw std::runtime_error(std::string("gsl error: ")
            + gsl_strerror(status));
      }
      return result;
   }


   double flrw_nat_A2::dl(double z) const
   {
    const double arg0 = acos(1 - m_crit) / 3.;
    const double arg1 = m_om / abs(m_ok);
    const double y1 = (-1. + cos(arg0) + M_SQRT3 * sin(arg0)) / 3.;
    const double y2 = (-1. - 2. * cos(arg0)) / 3.;
    const double y3 = (-1. + cos(arg0) - M_SQRT3 * sin(arg0)) / 3.;
    const double g = 2. / sqrt(y1 - y2);
    const double k = sqrt((y1 - y3) / (y1 - y2));
    const double phi = asin(sqrt((y1 - y2) / ((1 + z) * arg1 + y1)));
    const double phi0 = asin(sqrt((y1 - y2) / (arg1 + y1)));
    return (1. + z) / m_sqok * sin(g * (ellint_1(k, phi0) - ellint_1(k,
        phi)));
   }

   double flrw_nat_A2_1::age(double z) const
   {
     const double arg = (1 + z) * m_om / m_ok;
     return -(M_SQRT3 * log(1 - 2 / (sqrt(1. / 3. - arg) + 1)) + log(1 + 2
                    / (sqrt(1 - 3 * arg) - 1))) / sqrt(m_ov);
   }

   double flrw_nat_A2_2::age(double z) const
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

   // CASE A1
   double flrw_nat_A1::age(double z) const
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


  } // namespace impl

} // namespace milia

