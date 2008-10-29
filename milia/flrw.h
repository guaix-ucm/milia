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

#ifndef MILIA_FLRW_H
#define MILIA_FLRW_H

namespace milia
{
    namespace metrics
    {

        /**
         * The Friedmann-Lemaître-Robertson-Walker metric
         * 
         * This class computes the common cosmological distances and times.
         * It uses elliptical functions from gsl.
         * It is based on the paper <a href="http://xxx.unizar.es/abs/astro-ph/9905116">%astro-ph/9905116</a> 
         * for the relations between distances.
         * The age is computed from 
         * <a href="http://xxx.unizar.es/abs/astro-ph/0003463">%astro-ph/0003463</a> with elliptical functions.
         * Distances are computed from the luminosity distance using 
         * <a href="http://xxx.unizar.es/abs/astro-ph/0002334">%astro-ph/0002334</a> 
         * without inhomogeneities.
         */
        class flrw
        {
        public:

            /**
             * No Big Bang: \f[\Omega_v > 4 \Omega_m [f(\frac{1}{3}f^{-1}(\Omega_m^{-1} - 1))]^3 \f]
             * where \f[ f = \cos \Omega_m > 0.5 \cosh \Omega_m < 0.5 \f]
             * 
             * Recollapse: \f[\Omega_v < 0\f] or \f[\Omega_v > 0 \Omega_m > 1 
             * \Omega_v < 4 \Omega_m [\cos(\frac{1}{3}\cos^{-1}(\Omega_m^{-1} - 1) + \frac{4\pi}{3})]^3 \f]
             * 
             * From Cosmological Physics, Peacock pags 82-83          
             * 
             * @pre hubble_constant > 0 matter_density >= 0 lambda_density >= 0
             * @pre the parameters allow a Big-Bang to ocurr and the Universe doesn't recollapse
             * @throws milia::recollapse
             * @throws milia::no_big_bang
             * @throws milia::exception
             * 
             */
            flrw(double hubble, double matter, double vacuum);

            bool set_hubble(double hubble_parameter);

            inline double get_hubble() const
            {
                return m_hu;
            }

            inline double get_matter() const
            {
                return m_om;
            }

            bool set_matter(double matter_density);

            inline double get_vacuum() const
            {
                return m_ov;
            }

            bool set_vacuum(double vacuum_energy_density);

            double hubble(double z) const;

            double helper(double z) const;

            /**
             * Comoving distance (line of sight)
             * \f[
             * dc(z)=\frac{1}{H_0}\int_0^z \frac{dt}{\sqrt{\Omega_m(1+t)^3+\Omega_k(1+t)^2+\Omega_v}}
             * \f]
             */                  
            double dc(double z) const;

            /**
             * Comoving distance (transverse)
             */
            double dm(double z) const;

            /**
             * Angular distance
             */
            double da(double z) const;

            /**
             * Luminosity distance
             */                    
            double dl(double z) const;

            /**
             * Distance modulus
             */
            double DM(double z) const;

            /**
             * Comoving volume per solid angle
             */
            double vol(double z) const;

            /**
             * Current age of the Universe
             */
            double age() const;

            /**
             * Age of the Universe
             */
            double age(double z) const;

            /**
             * Look-back time
             */
            double lt(double z) const;

            double dc(double z, double dm) const;
            double dm(double z, double dl) const;
            double da(double z, double dl) const;
            double vol(double z, double dm) const;

        private:

            static const double ms_hubble_radius;

            static const double ms_hubble_time;

            /**
             * Hubble parameter
             */
            double m_hu;

            /**
             * Matter density
             */
            double m_om;

            /**
             * Vaccum energy density
             */
            double m_ov;

            /**
             * Separation parameter
             */
            double m_b;

            /**
             * Curvature parameter
             * 
             * m_om + m_ov + m_ok = 1
             */
            double m_ok;

            /**
             * Sign of the curvature parameter
             */
            int m_kap;

            double m_r_h;

            double m_t_h;

            double m_universe_age;

            enum CASES {
                OM_OV_0 = 1, //om = ov=0
                OV_1, //ov=0 0<om<1
                OV_2, //ov=0 om>1
                OV_3, //ov=0 om=1
                OM, //om=0
                A1, //om+ov!=1 b<0 || b>2
                A2_1, //om+ov!=1 b=2
                A2_2, //om+ov!=1 0<b<2
                OM_OV_1, //om+ol=1
                NO_CASE = 0,
            } m_case;

            CASES check() const;

            double tolz(double z) const;
            double tomz(double z) const;
            double ta1(double z) const;
            double ta2(double z) const;
            double tb(double z) const;
            double integrate_time(double) const;

            // Cache
            mutable struct
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
            }
            cache;

            void init_cache(double z) const;
            void scale_cache(double rh, double th) const;
            bool use_cache(double z) const;
        };

    } // namespace metrics

} // namespace milia


#endif /* MILIA_FLRW_H */
