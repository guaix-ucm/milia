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

#ifndef MILIA_FLRW_H
#define MILIA_FLRW_H

#include <string>
#include <ostream>

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
class flrw
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
     * Set the Hubble parameter \f[ H_0\f]
     *
     * @pre Hubble parameter > 0
     * @param hubble Hubble parameter in \f$ km\ s^{-1}\ Mpc^{-1} \f$
     * @return True if value is acceptable
     * @throw milia::exception
     */
    bool set_hubble(double hubble);

    /**
     * Get the value of the Hubble parameter in \f$ km\ s^{-1}\ Mpc^{-1} \f$.
     */
    inline double get_hubble() const
    {
        return m_hu;
    }

    /**
     * Get the value of the matter density \f[\Omega_m \f]
     */
    inline double get_matter() const
    {
        return m_om;
    }

    /**
     * Set the value of the matter density  \f[\Omega_m \f]
     *
     * @param M matter density
     * @return True if value acceptable
     * @throws milia::recollapse
     * @throws milia::no_big_bang
     * @throws milia::exception
     */
    bool set_matter(double M);

    /**
     * Get the value of the vacuum energy density \f[ \Omega_v \f]
     *
     */
    inline double get_vacuum() const
    {
        return m_ov;
    }

    /**
     * Set the value of the vacuum energy density \f[ \Omega_v \f]
     *
     * @param L vacuum energy density
     * @return True if value acceptable
     * @throws milia::recollapse
     * @throws milia::no_big_bang
     *
     */
    bool set_vacuum(double L);

    /**
     * Computes the Hubble parameter at redshift z
     *
     * @param z redshift
     * @return the Hubble parameter at the given redshift
     */
    double hubble(double z) const;

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
     * Comoving distance (line of sight) in Mpc
     * \f[
     * D_c(z)=\left\{
     * \begin{array}{rl}
     * D_m(z) & \textrm{if } \Omega_k = 0 \\
         * \frac{1}{\sqrt{|\Omega_k|}} \frac{c}{H_0} \sinh^{-1}(\sqrt{|\Omega_k|} \frac{H_0}{c} D_m(z) )   & \textrm{if } \Omega_k > 0 \\
         * \frac{1}{\sqrt{|\Omega_k|}} \frac{c}{H_0} \sin^{-1}(\sqrt{|\Omega_k|} \frac{H_0}{c} D_m(z))   & \textrm{if } \Omega_k < 0 \\
         * \end{array} \right.
     * \f]
     *
     * @param z redshift
     * @param dm comoving transverse distance in Mpc
     * @return line of sight comoving distance in Mpc
     *
     */
    double dc(double z, double dm) const;

    /**
     * Comoving distance (transverse) in Mpc
     *
     * @param z redshift
     * @return transverse comoving distance in Mpc
     *
     */
    double dm(double z) const;

    /**
     * Comoving distance (transverse) in Mpc
     *
     * \f[
     * D_m(z) = \frac{1}{1 + z} D_l(z)
     * \f]
     *
     * @param z redshift
     * @param z redshift
     * @param dl luminosity distance in Mpc
     * @return transverse comoving distance in Mpc
     *
     */
    double dm(double z, double dl) const;

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
     * Angular distance in Mpc
     * \f[
     * D_a(z) = \frac{1}{(1 + z)^2} D_l(z)
     * \f]
     *
     * @param z redshift
     * @param dl luminosity distance in Mpc
     * @return the angular distance in Mpc
     */
    double da(double z, double dl) const;

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
     * Comoving volume per solid angle
     *
     * @param z redshift
     * @param dm comoving transverse distance in Mpc
     * @return comoving volume in \f$ Mpc^3\f$ per solid angle
     */
    double vol(double z, double dm) const;

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
     * Transform angular sizes in arc seconds into parsecs
     *
     * @param z redshift
     * @param arcsec apparent size on the sky (in arcsec)
     * @return distance in parsec.
     *
     */
    double arcsec2pc(double z, double arcsec) const;

    /**
     * Transform physical lengths in parcsecs into angular sizes in arcseconds
     *
     * @param z redshift
     * @param pc distance in parsecs.
     * @return apparent size on the sky (in arcsec)
     */
    double pc2arcsec(double z, double pc) const;

private:
    static const double ms_hubble_radius = 2.99792e5; // Hubble Radius in Mpc
    static const double ms_hubble_time = 9.78e2; // Hubble time in Gyr

    /**
     * Hubble parameter
     */
    double m_hu;

    /**
     * Matter density
     */
    double m_om;

    /**
     * Vacuum energy density
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

    enum ComputationCases
    {
        NO_CASE, OM_OV_0, //om = ov=0
        OV_1, //ov=0 0<om<1
        OV_2, //ov=0 om>1
        OV_3, //ov=0 om=1
        OM, //om=0
        A1, //om+ov!=1 b<0 || b>2
        A2_1, //om+ov!=1 b=2
        A2_2, //om+ov!=1 0<b<2
        OM_OV_1,
    //om+ol=1
    };

    ComputationCases m_case;
    ComputationCases check() const;

    double tolz(double z) const;
    double tomz(double z) const;
    double ta1(double z) const;
    double ta2(double z) const;
    double tb(double z) const;
    double ti(double z) const;
};

class flrw_cache: public flrw
{
public:
    flrw_cache(double hubble, double matter, double vacuum);
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

std::ostream& operator<<(std::ostream& os, milia::metrics::flrw& iflrw);

#endif /* MILIA_FLRW_H */
