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

//#include <cmath>

#include <boost/math/special_functions/pow.hpp>

#include "flrw_cached.h"

using boost::math::pow;

namespace milia
{
  namespace metrics
  {

    flrw_cached::flrw_cached(double hubble, double matter, double vacuum) :
      metrics::flrw(hubble, matter, vacuum)
    {
    }

    bool flrw_cached::set_hubble(double hubble_parameter)
    {
      metrics::flrw::set_hubble(hubble_parameter);
      m_cache.recompute();
      return true;
    }

    bool flrw_cached::set_matter(double matter_density)
    {
      metrics::flrw::set_matter(matter_density);
      m_cache.recompute();
      return true;
    }

    bool flrw_cached::set_vacuum(double vacuum_energy_density)
    {
      metrics::flrw::set_vacuum(vacuum_energy_density);
      m_cache.recompute();
      return true;
    }

    double flrw_cached::hubble(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.hubble;
    }

    double flrw_cached::dc(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.dc;
    }

    double flrw_cached::da(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.da;
    }

    double flrw_cached::dm(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.dm;
    }

    double flrw_cached::dl(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.dl;
    }

    double flrw_cached::DM(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.DM;
    }

    double flrw_cached::vol(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.vol;
    }

    double flrw_cached::age() const
    {
      return m_cache.age_0;
    }

    double flrw_cached::age(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.age;
    }

    double flrw_cached::lt(double z) const
    {
      if (!m_cache.can_use(z))
      {
        m_cache.initialize(*this, z);
      }
      return m_cache.lt;
    }

    void flrw_cached::Cache::initialize(const metrics::flrw& metric, double zz)
    {
      // redshift
      z = zz;
      // hubble
      hubble = metric.hubble(z);
      // dl
      dl = metric.dl(z);
      // dm
      dm = metric.dm(z, dl);
      // da
      da = metric.da(z, dl);
      // dc
      dc = metric.dc(z, dm);
      // DM (check z > 0)
      if (z > 0)
        DM = 5. * log10(dl) + 25.;
      else
        DM = -1.;
      // Vol
      vol = metric.vol(z, dm);
      // Age
      age_0 = metric.age();
      age = metric.age(z);
      // Look-back time
      lt = age_0 - age;
    }

    bool flrw_cached::Cache::can_use(double zz) const
    {
      if (z == zz)
        return true;
      return false;
    }

    void flrw_cached::Cache::scale(double rh, double th)
    {
      dl /= rh;
      dm /= rh;
      da /= rh;
      dc /= rh;
      if (z > 0)
        DM = 5. * log10(dl) + 25.;
      else
        DM = -1.;
      vol /= pow<3> (rh);
      age /= th;
      lt /= th;
    }

    void flrw_cached::Cache::recompute()
    {
      z = -1;
    }

  } //namespace metrics

} //namespace milia
