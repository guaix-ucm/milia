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

#include <milia/flrw.h>

namespace milia
{
  namespace metrics
  {

    /**
     * The Friedmann-Lemaître-Robertson-Walker metric
     *
     * This class represents a FLRW metric. Its methods compute the
     * common cosmological distances and times.
     * It uses elliptical functions from boost.
     * It is based on the paper <a href="http://xxx.unizar.es/abs/astro-ph/9905116">%astro-ph/9905116</a>
     * for the relations between distances.
     * The age and distance luminosity are computed from
     * <a href="http://xxx.unizar.es/abs/astro-ph/0003463">%astro-ph/0003463</a> with elliptical functions.
     * Distances are computed from the luminosity distance using
     * <a href="http://xxx.unizar.es/abs/astro-ph/0002334">%astro-ph/0002334</a>
     * without inhomogeneities.
     */
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
