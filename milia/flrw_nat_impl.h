/*
 * Copyright 2008-2013 Sergio Pascual
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

#ifndef MILIA_FLRW_NAT_IMPL_H
#define MILIA_FLRW_NAT_IMPL_H

#include <cmath>

using std::abs;

namespace milia
{
  namespace impl
  {
    class flrw_nat_impl
    {
      public:
        flrw_nat_impl(double matter, double vacuum) : 
          m_om(matter), 
          m_ov(vacuum), 
          m_ok(1 - m_om - m_ov), 
          m_sqok(sqrt(abs(m_ok)))
        {
          m_kap = m_ok > 0 ? -1 : 1;
        }

        static flrw_nat_impl* construct(double matter, double vacuum);

        virtual double dc(double z) const = 0;

        double get_matter() const
        {
          return m_om;
        }

        double get_vacuum() const
        {
          return m_om;
        }
        double dm(double z) const 
        {
          return dl(z) / (1 + z);
        }
        
        double da(double z) const 
        {
          return dl(z) / ((1 + z) * (1 + z));
        }

        virtual double dl(double z) const = 0;

        virtual double vol(double z) const = 0;
        // double age() const = 0;
        virtual double age(double z) const = 0;
        //double lt(double z) const = 0;
      protected:
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

      private:
        // Current Universe age (may be infinity in certain models)
        double m_uage;

    };

  } // namespace impl

} // namespace milia

#endif /* MILIA_FLRW_NAT_IMPL_H */
