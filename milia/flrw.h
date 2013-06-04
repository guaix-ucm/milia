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

#ifndef MILIA_FLRW_H
#define MILIA_FLRW_H

#include <milia/flrw_nat.h>

#include <string>
#include <ostream>
#include <cmath>

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
    class flrw : public flrw_nat
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
         * @param hubble a float, it's the Hubble parameter in \f$ km\ s^{-1}\ Mpc^{-1} \f$
         * @param matter a float, it's the matter density (dimensionless)
         * @param vacuum a float, it's the vaccum energy density (dimensionless)
         * @throws milia::recollapse
         * @throws milia::no_big_bang
         * @throws milia::exception
         *
         */
        flrw(double hubble, double matter, double vacuum);

        /**
         * Set the Hubble parameter \f[ H_0\f]
         *
         * @pre Hubble parameter > 0
         * @param hubble Hubble parameter in \f$ km\ s^{-1}\ Mpc^{-1} \f$
         * @throw milia::exception
         */
        void set_hubble(double hubble);

        /**
         * Get the value of the Hubble parameter in \f$ km\ s^{-1}\ Mpc^{-1} \f$.
         */
        double get_hubble() const;

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
         * Distance modulus \f$ DM = 5 log(\frac{D_l}{10\ pc}) \f$
         *
         * @param z redshift
         * @return distance modulus in mag
         *
         */
        double DM(double z) const;

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
         * String with caracteristics of the FLRW universe (Hubble parameter,
         * Matter density, vacuum energy density.
         *
         * @return a string
         */
        std::string to_string() const;

        /**
         * Factor to transform angular sizes in arc seconds into parsecs
         *
         * @param z redshift
         * @return scale factor.
         *
         */
        double angular_scale(double z) const;

      private:
        // Hubble Radius in Mpc for H = 1 km s^-1
        static const double ms_hubble_radius = 299792.458;
        // Hubble time in Gyr for H = 1 km s^-1
        static const double ms_hubble_time = 977.792222;

        // Hubble parameter
        double m_hu;
        // Current Hubble radius
        double m_r_h;
        // Current Hubble time
        double m_t_h;
    };

    inline double flrw::get_hubble() const
    {
      return m_hu;
    }

    inline double flrw::get_hubble(double z) const
    {
      return m_hu * flrw_nat::get_hubble(z);
    }

    inline double flrw::lt(double z) const
    {
      return m_t_h * flrw_nat::lt(z);
    }

    inline double flrw::dc(double z) const
    {
      return m_r_h * flrw_nat::dc(z);
    }

    inline double flrw::dm(double z) const
    {
      return m_r_h * flrw_nat::dm(z);
    }

    inline double flrw::da(double z) const
    {
      return m_r_h * flrw_nat::da(z);
    }

    inline double flrw::DM(double z) const
    {
      return 5 * log10(dl(z)) + 25;
    }

    inline double flrw::vol(double z) const
    {
      return m_r_h * m_r_h * m_r_h * flrw_nat::vol(z);
    }

    inline double flrw::dl(double z) const
    {
      return m_r_h * flrw_nat::dl(z);
    }

    inline double flrw::age() const
    {
      return m_t_h * flrw_nat::age();
    }

    inline double flrw::age(double z) const
    {
      return m_t_h * flrw_nat::age(z);
    }

} // namespace milia

std::ostream& operator<<(std::ostream& os, milia::flrw& iflrw);

#include <milia/flrw_nat_impl.h>


namespace milia
{
 namespace rei {
    class flrw_nat
    {
      public:
        flrw_nat(double matter, double vacuum);
        static bool does_recollapse(double matter, double vacuum);

        double get_matter() const;
        double get_vacuum() const;
        double get_hubble(double z) const;

        double dc(double z) const;
        double dm(double z) const;
        double da(double z) const;
        double dl(double z) const;
        double vol(double z) const;
        double age(double z) const;
        double lt(double z) const;
    const char* model() const
    {
      return m_impl->model();
    }
      private:
        std::auto_ptr<impl::flrw_nat_impl> m_impl;
    };

    inline double flrw_nat::get_matter() const
    {
      return m_impl->get_matter();
    }

    inline double flrw_nat::get_vacuum() const
    {
      return m_impl->get_vacuum();
    }

    inline double flrw_nat::dc(double z) const
    {
      return m_impl->dc(z);
    }

    inline double flrw_nat::dl(double z) const
    {
      return m_impl->dl(z);
    }

    inline double flrw_nat::dm(double z) const
    {
      return m_impl->dm(z);
    }

    inline double flrw_nat::da(double z) const
    {
      return m_impl->da(z);
    }

    inline double flrw_nat::vol(double z) const
    {
      return m_impl->vol(z);
    }

    inline double flrw_nat::age(double z) const
    {
      return m_impl->age(z);
    }


    class flrw
    {
      public:
        flrw(double hubble, double matter, double vacuum) :
          m_nat(matter, vacuum),
          m_hu(hubble),
          m_r_h(ms_hubble_radius / m_hu),
          m_t_h(ms_hubble_time / m_hu)
      {}
 
      double da(double z) const
      {
        return m_r_h * m_nat.da(z);
      }

      double dl(double z) const
      {
        return m_r_h * m_nat.dl(z);
      }

      double dc(double z) const
      {
        return m_r_h * m_nat.dc(z);
      }

      double dm(double z) const
      {
        return m_r_h * m_nat.dm(z);
      }

      double vol(double z) const
      {
        return m_r_h * m_r_h * m_r_h * m_nat.vol(z);
      }

      double age(double z) const
      {
        return m_t_h * m_nat.age(z);
      }
    const char* model() const
    {
      return m_nat.model();
    }
      private:
        flrw_nat m_nat;
        static const double ms_hubble_radius = 299792.458;
        static const double ms_hubble_time = 977.792222;
        double m_hu;
        double m_r_h;
        double m_t_h;
    };

  } // namespace rei
} // namespace milia

#endif /* MILIA_FLRW_H */
