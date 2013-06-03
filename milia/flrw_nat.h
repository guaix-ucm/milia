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

#ifndef MILIA_FLRW_NAT_H
#define MILIA_FLRW_NAT_H

#include <string>
#include <ostream>
#include <memory>

#include "flrw_nat_impl.h"

namespace milia
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
    class flrw_nat
    {
      public:

        /**
         * No Big Bang: \f[\Omega_v > 4 \Omega_m [f(\frac{1}{3}f^{-1}(\Omega_m^{-1} - 1))]^3 \f]
         * where \f[ f = \cos\ \textrm{if}\ \Omega_m > 0.5\ \textrm{and}\ f = \cosh\ \textrm{if}\ \Omega_m < 0.5 \f]
         *
         * Recollapse: \f[\Omega_v < 0\f] or \f[\Omega_v > 0,\ \Omega_m > 1\ \textrm{and}\
         * \Omega_v < 4 \Omega_m [\cos(\frac{1}{3}\cos^{-1}(\Omega_m^{-1} - 1) + \frac{4\pi}{3})]^3 \f]
         *
         * From Cosmological Physics, Peacock pags 82-83
         *
         * @pre Hubble parameter > 0 matter_density >= 0 lambda_density >= 0
         * @pre the parameters allow a Big-Bang to occur and the Universe doesn't recollapse
         * @param matter a float, it's the matter density (dimensionless)
         * @param vacuum a float, it's the vaccum energy density (dimensionless)
         * @throws milia::recollapse
         * @throws milia::no_big_bang
         * @throws milia::exception
         *
         */
        flrw_nat(double matter, double vacuum);

        /**
         * Checks whether the Universe recollapses or not
         * with the given parameters.
         *
         * Recollapse occurs if : \f[\Omega_v < 0\f] or \f[\Omega_v > 0,\ \Omega_m > 1\ \textrm{and}\
         * \Omega_v < 4 \Omega_m [\cos(\frac{1}{3}\cos^{-1}(\Omega_m^{-1} - 1) + \frac{4\pi}{3})]^3 \f]
         *
         * From Cosmological Physics, Peacock pags 82-83
         *
         * @param matter Matter density
         * @param vacuum Vacuum density
         * @return True if the Universe recollapses, false otherwise
         */
        static bool does_recollapse(double matter, double vacuum);

        /**
         * Get the value of the matter density \f[\Omega_m \f]
         */
        double get_matter() const;

        /**
         * Set the value of the matter density  \f[\Omega_m \f]
         *
         * @param matter matter density
         * @throws milia::recollapse
         * @throws milia::no_big_bang
         * @throws milia::exception
         */
        void set_matter(double matter);

        /**
         * Get the value of the vacuum energy density \f[ \Omega_v \f]
         *
         */
        double get_vacuum() const;

        /**
         * Set the value of the vacuum energy density \f[ \Omega_v \f]
         *
         * @param vacuum vacuum energy density
         * @throws milia::recollapse
         * @throws milia::no_big_bang
         *
         */
        void set_vacuum(double vacuum);

        /**
         * Computes the Hubble parameter at redshift z
         *
         * @param z redshift
         * @return the Hubble parameter at the given redshift
         */
        double get_hubble(double z) const;

        /**
         * Comoving distance (line of sight) in Mpc
         * \f[
         * D_c(z)=\frac{c}{H_0}\int_0^z \frac{dt}{\sqrt{\Omega_m(1+t)^3+\Omega_k(1+t)^2+\Omega_v}}
         * \f]
         *
         * @param z redshift
         * @return line of sight comoving distance in Mpc
         *
         */
        double dc(double z) const;

        /**
         * Comoving distance (transverse) in Mpc
         *
         * @param z redshift
         * @return transverse comoving distance in Mpc
         *
         */
        double dm(double z) const;

        /**
         * Angular distance in Mpc
         * \f[
         * D_a(z) = \frac{1}{1 + z} D_m(z)
         * \f]
         *
         * @param z redshift
         * @return the angular distance in Mpc
         */
        double da(double z) const;

        /**
         * Luminosity distance in Mpc
         *
         * @param z redshift
         * @return the luminosity distance in Mpc
         *
         */
        double dl(double z) const;

        /**
         * Comoving volume per solid angle
         *
         * @param z redshift
         * @return comoving volume in \f$ Mpc^3\f$ per solid angle
         */
        double vol(double z) const;

        /**
         * Current age of the Universe
         */
        double age() const;

        /**
         * Age of the Universe in Gyr
         *
         * @param z redshift
         * @return Age of the Universe in Gyr
         *
         */
        double age(double z) const;

        /**
         * Look-back time in Gyr
         *
         * @param z redshift
         * @return Look-back time in Gyr
         */
        double lt(double z) const;

        /**
         * String with characteristics of the FLRW universe (Hubble parameter,
         * Matter density, vacuum energy density.
         *
         * @return a string
         */
        std::string to_string() const;

      private:

        // Matter density
        double m_om;

        // Vacuum energy density
        double m_ov;

        // Critical parameter
        double m_crit;

        // Curvature parameter
        // m_om + m_ov + m_ok = 1
        double m_ok;
        // Square root of abs(m_ok)
        double m_sqok;
        // Negative of the sign of the curvature parameter
        short m_kap;
        // Current Universe age (may be infinity in certain models)
        double m_uage;

        enum ComputationCases
        {
          NO_CASE, // error condition
          OM_OV_0, //om = ov = 0
          OV_1, //ov = 0 0 < om < 1
          OV_2, //ov = 0 om > 1
          OV_EDS, //ov = 0 om = 1, Einstein-de Sitter Universe
          OM, //om = 0 0 < ov < 1
          OM_DS, //om = 0 ov = 1, de Sitter Universe
          OM_OV_1, //om + ol = 1
          A1, //om+ov != 1 b < 0 || b > 2
          A2_1, //om+ov != 1 b = 2
          A2_2, //om+ov != 1 0 < b < 2
        };

        ComputationCases m_case;
        ComputationCases select_case() const;

        // Distances and volumes from other distance
        double da(double z, double dl) const;
        double dc(double z, double dm) const;
        double dm(double z, double dl) const;
        double vol(double z, double dm) const;

        // Methods to compute the age in different cases
        double tolz(double z) const;
        double tomz(double z) const;
        double ta1(double z) const;
        double ta2(double z) const;
        double tb(double z) const;
        double ti(double z) const;
    };

    inline double flrw_nat::get_matter() const
    {
      return m_om;
    }

    inline double flrw_nat::get_vacuum() const
    {
      return m_ov;
    }

    inline double flrw_nat::dc(double z) const
    {
      return dc(z, dm(z));
    }

    inline double flrw_nat::dm(double z) const
    {
      return dm(z, dl(z));
    }

    inline double flrw_nat::dm(double z, double dl) const
    {
      return dl / (1. + z);
    }

    inline double flrw_nat::da(double z) const
    {
      return da(z, dl(z));
    }

    inline double flrw_nat::da(double z, double dl) const
    {
      return dm(z, dl) / (1 + z);
    }

    inline double flrw_nat::vol(double z) const
    {
      return vol(z, dm(z, dl(z)));
    }

} // namespace milia

std::ostream& operator<<(std::ostream& os, milia::flrw_nat& iflrw);

namespace milia
{
    class flrw_nat_new
    {
      public:
        flrw_nat_new(double matter, double vacuum);
        static bool does_recollapse(double matter, double vacuum);

        double get_matter() const;
        double get_vacuum() const;
        double get_hubble(double z) const;

        double dc(double z) const;
        double dm(double z) const;
        double da(double z) const;
        double dl(double z) const;
        double vol(double z) const;
        double age() const;
        double age(double z) const;
        double lt(double z) const;
      private:
        std::auto_ptr<impl::flrw_nat_impl> m_impl;

    };

    inline double flrw_nat_new::get_matter() const
    {
//      return m_impl->get_matter();
      return 0.0;
    }

    inline double flrw_nat_new::get_vacuum() const
    {
      return 0.0;
      // return m_impl->get_vacuum();
    }

    inline double flrw_nat_new::dc(double z) const
    {
      return m_impl->dc(z);
    }

    inline double flrw_nat_new::dl(double z) const
    {
      return m_impl->dl(z);
    }

    inline double flrw_nat_new::dm(double z) const
    {
      return m_impl->dm(z);
    }

    inline double flrw_nat_new::da(double z) const
    {
      return m_impl->da(z);
    }

    inline double flrw_nat_new::vol(double z) const
    {
      return m_impl->vol(z);
    }

} // namespace milia

#endif /* MILIA_FLRW_NAT_H */
