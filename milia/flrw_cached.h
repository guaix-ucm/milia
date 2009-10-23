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

#ifndef MILIA_FLRW_CACHED_H
#define MILIA_FLRW_CACHED_H

#include <string>
#include <ostream>

#include <milia/flrw.h>

namespace milia
{
  namespace metrics
  {
    class flrw_cached: public flrw
    {
      public:
        flrw_cached(double hubble, double matter, double vacuum);
        bool set_hubble(double hubble_parameter);
        bool set_matter(double matter_density);
        bool set_vacuum(double vacuum_energy_density);
        double hubble(double z) const;

        double dc(double z) const;
        double dm(double z) const;
        double da(double z) const;
        double dl(double z) const;
        double DM(double z) const;
        double vol(double z) const;
        double age() const;
        double age(double z) const;
        double lt(double z) const;

      private:

        struct Cache
        {
            double z;
            double hubble;
            double dc;
            double dm;
            double da;
            double dl;
            double DM;
            double vol;
            double lt;
            double age;
            double age_0;
            void recompute();
            bool can_use(double z) const;
            void scale(double rh, double th);
            void initialize(const metrics::flrw& metric, double z);
        };

        mutable Cache m_cache;
    };

  } // namespace metrics

} // namespace milia

#endif /* MILIA_FLRW_CACHED_H */
