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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "flrw.h"
#include "exception.h"

#include <cmath>
#include <cfloat>
#include <iostream>

#define M_SQRT3 1.73205080757
#define PREC GSL_PREC_SINGLE
#define EPS FLT_EPSILON

#include <gsl/gsl_sf_ellint.h>

namespace milia {

namespace metrics {

double flrw::distance_luminosity(double z) const {
	switch (m_case) {
	case OM_OL_0:
		return 0.5*z*(z+2);
		break;
	case OL_1:
	case OL_2:
	case OL_3:
		return m_r_h*2.*((2.-m_om*(1.-z)-(2.-m_om)*sqrt(1+m_om*z)))/(m_om*m_om);
		break;
	case OM: {
		//	cerr<<"No implementado todavia\n";
		return 0;
	}
		break;
	case A1:
		//om+ol!=1 b<0 || b>2
	{
		const double v=pow(m_kap*(m_b-1.)+sqrt(m_b*(m_b-2.)), 1./3.);
		const double y=(-1.+m_kap*(v+1./v))/3.;
		const double A=sqrt(y*(3.*y+2.));
		const double g=1./sqrt(A);
		const double k=sqrt(0.5+0.25*g*g*(v+1./v));
		const double sup=m_om/fabs(m_ok);
		const double phi=acos(((1.+z)*sup+m_kap*y-A)/((1.+z)*sup+m_kap*y+A));
		const double phi0=acos((sup+m_kap*y-A)/(sup+m_kap*y+A));
		if (m_b>0)
			return m_r_h*(1+z)/sqrt(fabs(m_ok))*sin(g*(gsl_sf_ellint_F(phi0, k,
					PREC)-gsl_sf_ellint_F(phi, k, PREC)));
		else
			return m_r_h*(1+z)/sqrt(fabs(m_ok))*sinh(g*(gsl_sf_ellint_F(phi0,
					k, PREC)-gsl_sf_ellint_F(phi, k, PREC)));
	}
	case A2_1:
		// b=2
	case A2_2:
		// 0<b<2
	{
		const double arg0=acos(1.-m_b)/3.;
		const double arg1=m_om/fabs(m_ok);
		const double y1=(-1.+cos(arg0)+M_SQRT3*sin(arg0))/3.;
		const double y2=(-1.-2.*cos(arg0))/3.;
		const double y3=(-1.+cos(arg0)-M_SQRT3*sin(arg0))/3.;
		const double g=2./sqrt(y1-y2);
		const double k=sqrt((y1-y3)/(y1-y2));
		const double phi=asin(sqrt((y1-y2)/((1+z)*arg1+y1)));
		const double phi0=asin(sqrt((y1-y2)/(arg1+y1)));
		return m_r_h*(1.+z)/sqrt(fabs(m_ok))*sin(g*(gsl_sf_ellint_F(phi0, k,
				PREC)-gsl_sf_ellint_F(phi, k, PREC)));
	}
	case OM_OL_1: {
		// m_om+ol=1
		const double k=0.9659258263;
		const double c1=1.316074013; //Raiz cuarta de 3
		const double arg0=pow((1/m_om-1), 1./3);
		const double down=1+(1+M_SQRT3)*arg0;
		const double up=1+(1-M_SQRT3)*arg0;
		const double phi=acos((z+up)/(z+down));
		const double phi0=acos(up/down);
		return m_r_h*(1+z)/(c1*sqrt(m_om)*sqrt(arg0))*(gsl_sf_ellint_F(phi0, k,
				PREC)-gsl_sf_ellint_F(phi, k, PREC));
	}
		break;
	}
	return -1;
}

} // namespace metrics

}
// namespace milia