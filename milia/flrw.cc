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

#include "flrw.h"
#include "exception.h"
#include "flrw_prec.h"

#include <cmath>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#ifndef HAVE_ASINH
#define asinh gsl_asinh
#endif

#ifndef HAVE_ATANH
#define atanh gsl_atanh
#endif

namespace milia
{
namespace metrics
{
using std::abs;

flrw::flrw(double h, double m, double v) :
    m_hu(h), m_om(m), m_ov(v), m_ok(1 - m_om - m_ov)
{
    gsl_set_error_handler_off();

    if (m_hu <= 0)
        throw milia::exception("Hubble constant <= 0 not allowed");

    // om < 0 not allowed
    if (m_om < -EPS)
        throw milia::exception("Matter density < 0 not allowed");

    //ol < 0 makes the universe recollapse
    if (m_ov < -EPS)
        throw milia::recollapse("The Universe recollapses"); // Recollapse

    m_b = -13.5 * gsl_pow_2(m_om) * m_ov / (gsl_pow_3(m_ok));

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

std::string flrw::to_string() const
{
    std::stringstream out;
    out << "flrw(hubble=" << m_hu << ", matter=" << m_om << ", vacuum=" << m_ov
            << ")";
    return out.str();
}

bool flrw::does_recollapse(double matter, double vacuum)
{
    if (vacuum < 0)
        return true;
    if (matter < 1)
        return false;
    const double critical = 4 * matter * gsl_pow_3(cos(1. / 3. * acos(1.
            / matter - 1.) + 4 * M_PI / 3.));
    if (vacuum > critical)
        return false;
    return true;
}

bool flrw::set_hubble(double H)
{
    if (H > 0)
    {
        m_hu = H;
        m_r_h = ms_hubble_radius / m_hu;
        m_t_h = ms_hubble_time / m_hu;
        m_universe_age = age();
        return true;
    }
    else
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

    return true;
}

bool flrw::set_vacuum(double L)
{
    if (L < -EPS)
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

    return true;
}

flrw::ComputationCases flrw::check() const
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
    if (l4 && m_om > 1 + EPS)
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
    if ((m_b <= EPS) || (m_b > 2 + EPS))
        return A1;
    // a2 b=2
    if ((m_b >= 2 - EPS) && (m_b <= 2 + EPS))
        return A2_1;
    // a2 0<b<2
    if ((m_b > EPS) && (m_b < 2 - EPS))
        return A2_2;
    //
    return NO_CASE;
}

double flrw::hubble(double z) const
{
    return m_hu
            * sqrt(m_om * gsl_pow_3(1 + z) + m_ok * gsl_pow_2(1 + z) + m_ov);
}

double flrw::lt(double z) const
{
    return m_universe_age - age(z);
}

double flrw::dc(double z) const
{
    return dc(z, dm(z));
}

double flrw::dc(double z, double dm) const
{
    if (m_ok == 0)
        return dm;
    else
    {
        double sqok = sqrt(abs(m_ok));
        if (m_ok > 0)
            return m_r_h / sqok * asinh(dm * sqok / m_r_h);
        else
            return m_r_h / sqok * asin(dm * sqok / m_r_h);
    }
}

double flrw::dm(double z) const
{
    return dm(z, dl(z));
}

double flrw::dm(double z, double dl) const
{
    return dl / (1. + z);
}

double flrw::da(double z) const
{
    return da(z, dl(z));
}

double flrw::da(double z, double dl) const
{
    return dl / gsl_pow_2(1 + z);
}

double flrw::DM(double z) const
{
    return 5 * log10(dl(z)) + 25;
}

double flrw::vol(double z) const
{
    return vol(z, dm(z, dl(z)));
}

double flrw::vol(double z, double dm) const
{
    if (abs(m_ok) < EPS)
    {
        return gsl_pow_3(dm) / 3.0;
    }
    else
    {
        const double r = m_hu / ms_hubble_radius * dm;
        const double sqrtok = sqrt(abs(m_ok));
        if (m_ok > EPS)
            return 0.5 * gsl_pow_3(m_r_h) * (r * sqrt(1 + m_om * r) - asinh(
                    sqrtok * r)) / sqrtok / m_ok;
        else
            return 0.5 * gsl_pow_3(m_r_h) * (asin(sqrtok * r) - r * sqrt(1
                    + m_om * r)) / sqrtok / m_ok;
    }
}

double flrw::arcsec2pc(double z, double arcsec) const
{
    // 206264.8062 converts arcscondss to radians
    return dl(z) * 1e6 / ((1 + z) * (1 + z)) / 206264.8062 * arcsec;
}

double flrw::pc2arcsec(double z, double pc) const
{
    // 206264.8062 converts arcscondss to radians
    return ((1 + z) * (1 + z)) / (dl(z) * 1e6) * 206264.8062 * pc;
}

flrw_cache::flrw_cache(double hubble, double matter, double vacuum) :
    metrics::flrw(hubble, matter, vacuum)
{
}

bool flrw_cache::set_hubble(double hubble_parameter)
{
    metrics::flrw::set_hubble(hubble_parameter);
    m_cache.recompute();
    return true;
}

bool flrw_cache::set_matter(double matter_density)
{
    metrics::flrw::set_matter(matter_density);
    m_cache.recompute();
    return true;
}

bool flrw_cache::set_vacuum(double vacuum_energy_density)
{
    metrics::flrw::set_vacuum(vacuum_energy_density);
    m_cache.recompute();
    return true;
}

double flrw_cache::hubble(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.hubble;
}

double flrw_cache::dc(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.dc;
}

double flrw_cache::da(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.da;
}

double flrw_cache::dm(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.dm;
}

double flrw_cache::dl(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.dl;
}

double flrw_cache::DM(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.DM;
}

double flrw_cache::vol(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.vol;
}

double flrw_cache::age() const
{
    return m_cache.age_0;
}

double flrw_cache::age(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.age;
}

double flrw_cache::lt(double z) const
{
    if (!m_cache.can_use(z))
    {
        m_cache.initialize(*this, z);
    }
    return m_cache.lt;
}

void flrw_cache::Cache::initialize(const metrics::flrw& metric, double zz)
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

bool flrw_cache::Cache::can_use(double zz) const
{
    if (z == zz)
        return true;
    return false;
}

void flrw_cache::Cache::scale(double rh, double th)
{
    dl /= rh;
    dm /= rh;
    da /= rh;
    dc /= rh;
    if (z > 0)
        DM = 5 * log10(dl) + 25;
    else
        DM = -1.;
    vol /= gsl_pow_3(rh);
    age /= th;
    lt /= th;
}

void flrw_cache::Cache::recompute()
{
    z = -1;
}

} //namespace metrics
} //namespace milia

std::ostream& operator<<(std::ostream& os, milia::metrics::flrw& iflrw)
{
    os << iflrw.to_string();
    return os;
}
