/*
 * Copyright 2008 Sergio Pascual
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

#include "flrw.h"
#include "exception.h"
#include <cmath>
#include <gsl/gsl_math.h>

#ifndef HAVE_ASINH
#define asinh gsl_asinh
#endif

#ifndef HAVE_ATANH
#define atanh gsl_atanh
#endif

#define PREC GSL_PREC_DOUBLE
#define EPS 1e-5

namespace milia
{
  namespace metrics
  {
    using std::abs;

    const double flrw::ms_hubble_radius=2.99792e5; // Hubble Radius in Mpc
    const double flrw::ms_hubble_time=9.78e2; // Hubble time in Gyr

    flrw::flrw(double h, double m, double v) :
      m_hu(h), m_om(m), m_ov(v), m_ok(1 - m_om - m_ov)
    {

      if (m_hu <= 0)
        throw milia::exception("Hubble constant <= 0 not allowed");

      //om < 0 not allowed
      if (m_om < -EPS)
        throw milia::exception("Matter density < 0 not allowed");

      //ol < 0 makes the universe recollapse
      if (m_ov < -EPS)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      m_b= -13.5 * gsl_pow_2(m_om) * m_ov / (gsl_pow_3(m_ok));

      if (m_om >= 1 && m_b <= 2)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      if (m_ov >= 1 && m_b <= 2)
        throw milia::no_big_bang("No Big Bang"); // No Big bang with this parameters

      m_kap = (m_ok > 0 ? -1 : 1);
      m_case = check();
      m_r_h = ms_hubble_radius / m_hu;
      m_t_h = ms_hubble_time / m_hu;
      m_universe_age = age();
    }

    bool flrw::does_recollapse(double matter, double vacuum)
    {
      if (vacuum < 0)
        return true;
      if (matter < 1)
        return false;
      const double critical = 4 * matter * gsl_pow_3(cos(1. / 3. * acos(1.
          / matter - 1.) + 4 * M_PI / 3.));
      if (vacuum> critical)
        return false;
      return true;
    }

    bool flrw::use_cache(double z) const
    {
      if (cache.z==z)
        return true;
      return false;
    }

    void flrw::init_cache(double z) const
    {
      // redshift
      cache.z=z;
      // hubble
      cache.hubble=m_hu*sqrt(m_om*(1+z)*(1+z)*(1+z)+m_ok*(1+z)*(1+z)+m_ov);
      // dl
      cache.dl = dl(z);
      // dm
      cache.dm = dm(z, cache.dl);
      // da
      cache.da = da(z, cache.dl);
      // dc
      cache.dc = dc(z, cache.dm);
      // DM (check z > 0)
      if (cache.z> 0)
        cache.DM = 5. * log10(cache.dl) + 25.;
      else
        cache.DM = -1.;
      // Vol
      cache.vol = vol(z, cache.dm);
      // Age
      cache.age = age(z);
      // Look-back time
      cache.lt = m_universe_age - cache.age;
    }

    void flrw::scale_cache(double rh, double th) const
    {
      cache.dl/=rh;
      cache.dm/=rh;
      cache.da/=rh;
      cache.dc/=rh;
      if (cache.z>0)
        cache.DM=5.*log10(cache.dl)+25.;
      else
        cache.DM=-1.;
      cache.vol/=gsl_pow_3(rh);
      cache.age/=th;
      cache.lt/=th;
    }

    bool flrw::set_hubble(double H)
    {
      if (H > 0)
      {
        m_hu = H;
        m_r_h = ms_hubble_radius / m_hu;
        m_t_h = ms_hubble_time / m_hu;
        m_universe_age = age();
        // Recomputes cache the next call
        cache.z=-1;
        return true;
      } else
      {
        throw milia::exception("Hubble constant <= 0 not allowed");
      }
    }

    bool flrw::set_matter(double M)
    {
      if (M < -EPS)
      {
        throw milia::exception("Matter density < 0 not allowed");
      }

      const double OK = 1 - m_ov - M;
      double B = -13.5 * gsl_pow_2(M) * m_ov / gsl_pow_3(OK);

      if (M >= 1 && B <= 2)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      if (m_ov >= 1 && B <= 2)
        throw milia::no_big_bang("No Big Bang"); // No Big bang with this parameters

      // Escalas de tamanyo y tiempo;
      m_om = M;
      m_ok = OK;
      m_b = B;
      m_kap = (m_ok > 0 ? -1 : 1);
      m_case = check();
      m_universe_age = age();

      // recompute cache in the next call
      cache.z = -1;
      return true;
    }

    bool flrw::set_vacuum(double L)
    {
      if (L<-EPS)
        throw milia::recollapse("The Universe recollapses"); // Recollapse        

      const double OK = 1 - m_om - L;
      double B = -13.5 * gsl_pow_2(m_om) * L / gsl_pow_3(OK);

      if (m_om >= 1 && B <= 2)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

      if (L >= 1 && B <= 2)
        throw milia::no_big_bang("No Big Bang"); // No Big bang with this parameters

      m_ov = L;
      m_ok = OK;
      m_b = B;
      m_kap = (m_ok > 0 ? -1 : 1);
      m_case = check();
      m_universe_age = age();

      // Recompute cache next time
      cache.z=-1;
      return true;
    }

    flrw::CASES flrw::check() const
    {
      const bool l3 = (abs(m_om) <= EPS);
      const bool l4 = (abs(m_ov) <= EPS);
      // om=ov=0
      if (l3 && l4)
        return OM_OV_0;
      // ov=0 0<om<1
      if (l4 && m_om < 1 - EPS)
        return OV_1;
      // ov=0 om>1
      if (l4 && m_om> 1 + EPS)
        return OV_2;
      // ov=0 om==1
      if (l4)
        return OV_3;
      // om=0
      if (l3)
        return OM;
      // om+ov==1
      if (abs(m_ok) <= EPS)
        return OM_OV_1;
      // om+ov!=1
      // a1 b<0 || b>2
      if ((m_b <= EPS) || (m_b> 2 + EPS))
        return A1;
      // a2 b=2
      if ((m_b >= 2 - EPS) && (m_b <= 2 + EPS))
        return A2_1;
      // a2 0<b<2
      if ((m_b> EPS) && (m_b < 2 - EPS))
        return A2_2;
      //
      return NO_CASE;
    }

    double flrw::hubble(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.hubble;
    }

    double flrw::helper(double z) const
    {
      return sqrt(m_om*gsl_pow_3(1+z)+m_ok*gsl_pow_2(1+z)+m_ov);
      // Which one is better?
      // return sqrt(gsl_pow_2(1.+z)*(1.+om*z)-z*ol*(2.+z))));
    }

    double flrw::lt(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.lt;
    }

    double flrw::dc(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.dc;
    }

    double flrw::dc(double z, double dm) const
    {
      if (abs(m_ok) < EPS)
        return dm;
      else
      {
        double sqok = sqrt(abs(m_ok));
        if (m_ok> EPS)
          return m_r_h / sqok * asinh(dm* sqok / m_r_h);
        else
          return m_r_h / sqok * asin(dm * sqok / m_r_h);
      }
    }

    double flrw::dm(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.dm;
    }

    double flrw::dm(double z, double dl) const
    {
      return dl / (1. + z);
    }

    double flrw::da(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.da;
    }

    double flrw::da(double z, double dl) const
    {
      return dl / gsl_pow_2(1. + z);
    }

    double flrw::DM(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.DM;
    }

    // volumen por estereoradianes
    double flrw::vol(double z) const
    {
      if (!use_cache(z))
      {
        init_cache(z);
      }
      return cache.vol;
    }

    double flrw::vol(double z, double dm) const
    {
      if (abs(m_ok) < EPS)
      {
        return gsl_pow_3(dm) / 3.0;
      } else
      {
        const double r = m_hu / ms_hubble_radius * dm;
        const double sqrtok = sqrt(abs(m_ok));
        if (m_ok> EPS)
          return 0.5 * gsl_pow_3(m_r_h) * (r * sqrt(1 + m_om * r) - asinh(sqrtok * r)) / sqrtok / m_ok;
        else
          return 0.5 * gsl_pow_3(m_r_h) * (asin(sqrtok*r) - r * sqrt(1 + m_om
              * r)) / sqrtok / m_ok;
      }
    }
  } //namespace metrics
} //namespace milia