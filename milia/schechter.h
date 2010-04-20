/*
 * Copyright 2008-2010 Sergio Pascual
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

#ifndef MILIA_SCHECHTER_H
#define MILIA_SCHECHTER_H

#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <ostream>

namespace milia
{
  namespace luminosity_functions
  {
    /**
     * A Schechter luminosity function
     */
    class schechter
    {
      public:

        /**
         * Template constructor with evolution
         */
        template<typename EvolP, typename EvolL, typename EvolA> schechter(
            const EvolP& phi_star, const EvolL& lum_star, const EvolA& alpha,
            double z) :
          m_phi_star(phi_star), m_lum_star(lum_star), m_alpha(alpha),
              m_current_z(z)
        {
        }

        /**
         * Constructor without evolution
         */
        schechter(double phi_star, double lum_star, double alpha);

        /**
         * Generic constructor with evolution
         */
        schechter(const boost::function<double(double)>& phi_star, const boost::function<double(double)>& lum_star, const boost::function<double(double)>& alpha, double z);

        /**
         * Constructor with evolution parametrized following Heyl et al 1997
         */
        schechter(double phi_star, double e_phi_star, double lum_star,
            double e_lum_star, double alpha, double e_alpha, double redshift);

        double function(double lum) const;
        double object_density(double lum1, double lum2) const;
        /*double luminosity_density(double lum1, double lum2) const;*/
        void evolve(double z);
        boost::tuple<double, double, double> parameters() const;
        std::string to_string() const;
      private:
        boost::function<double(double)> m_phi_star;
        boost::function<double(double)> m_lum_star;
        boost::function<double(double)> m_alpha;
        double m_current_z;
        double function_normalized(double x) const;
    };

  } // namespace luminosity_functions

} // namespace milia

std::ostream& operator<<(std::ostream& os,
		milia::luminosity_functions::schechter& ischech);

#endif /* MILIA_SCHECHTER_H */
