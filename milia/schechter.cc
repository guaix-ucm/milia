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

#include "schechter.h"

#include <cmath>
#include <sstream>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>



namespace milia
{
  namespace luminosity_functions
  {
    using namespace boost::lambda;
    typedef boost::function<double(double)> OOFun;

    schechter::schechter(double phi_star, double lum_star, double alpha) :
      m_phi_star((_1 = phi_star, _1)), m_lum_star((_1 = lum_star, _1)),
          m_alpha((_1 = alpha, _1)), m_current_z(0.)
    {
    }

    schechter::schechter(const OOFun& phi_star, const OOFun& lum_star,
        const OOFun& alpha, double z) :
      m_phi_star(phi_star), m_lum_star(lum_star), m_alpha(alpha),
          m_current_z(z)
    {
    }

    schechter::schechter(double phi_star, double e_phi_star, double lum_star,
        double e_lum_star, double alpha, double e_alpha, double z) :
      m_current_z(z)
    {
      const double p0 = phi_star / pow((1 + z), e_phi_star);
      const double l0 = lum_star / pow((1 + z), e_lum_star);

      m_phi_star = p0 * bind(pow, 1. + _1, e_phi_star);
      m_phi_star = l0 * bind(pow, 1. + _1, e_lum_star);
      m_alpha = alpha + e_alpha * (_1 - z);
    }

    std::string schechter::to_string() const
    {
    	std::stringstream out;
    	out << "schechter(" << "phi*=" << m_phi_star(m_current_z)
    	<< " ,lum*=" << m_lum_star(m_current_z)
    	<< " ,alpha=" << m_alpha(m_current_z)
    	<< " ,z=" << m_current_z << ")";
    	return out.str();
    }

    boost::tuple<double,double,double> schechter::parameters() const
    {
      return boost::make_tuple(m_phi_star(m_current_z),
          m_lum_star(m_current_z), m_alpha(m_current_z));
    }

    void schechter::evolve(double z)
    {
      m_current_z = z;
    }

    double schechter::function(double lum) const
    {
      const double x = lum / m_lum_star(m_current_z);
      const double alpha = m_alpha(m_current_z);
      const double phi = m_phi_star(m_current_z);
      return phi * pow(x, alpha) * std::exp(-x);
    }

    double schechter::function_normalized(double x) const
    {
      const double alpha = m_alpha(m_current_z);
      return std::pow<double,double>(x, alpha) * std::exp(-x);
    }

    double schechter::object_density(double lum1, double lum2) const
    {
      const double lum = m_lum_star(m_current_z);
      const double x1 = lum1 / lum;
      const double x2 = lum2 / lum;
      const double phi = m_phi_star(m_current_z);

      static const int NSIMPSON = 1001;
      const double h = (x2 - x1) / (NSIMPSON - 1);
      double s = 0;

      for (int i = 2; i < NSIMPSON; i++)
      {
        const int f = ((i % 2) == 0 ? 4 : 2);
        s += f * function_normalized((x1 + (i - 1) * h));
      }
      return phi * h * (function_normalized(x1) + function_normalized(x2) + s)
          / 3;

      /*return phi * (gsl_sf_gamma_inc(alpha + 1, x1) - gsl_sf_gamma_inc(alpha
       + 1, x2));*/
    }

  } // namespace luminosity_functions

}
// namespace milia

std::ostream& operator<<(std::ostream& os,
		milia::luminosity_functions::schechter& ischech)
{
	os << ischech.to_string();
	return os;
}
