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

#define M_SQRT3 1.73205080757
#define PREC GSL_PREC_SINGLE
#define EPS 1e-4

#include <cmath>
#include <cfloat>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>

#include "flrw.h"

/* 
 OJO. La funcion elliptica con tres argumento sigue distinto
 criterio en GSl y en NR
 ellpi(phi,n,k)=gsl_sf_ellint_P(phi,k,n)
 Ademas este pollo usa otro convenio con los signos, asi que:
 ellpi(phi,n,k)=gsl_sf_ellint_P(phi,k,-n)
 */

namespace {
double helper_fun_time(double z, void* pars) {
	milia::metrics::flrw* p_cosm = static_cast<milia::metrics::flrw*>(pars);
	return 1./((1+z)*(p_cosm->helper(z)));
}
}

namespace milia {
namespace metrics {

double flrw::age() const {
	switch (m_case) {
	case OM_OL_0:
		return m_t_h;
		break;
	case OL_1:
		return m_t_h*(1.0-m_om/sqrt(1.0-m_om)*atanh(sqrt(1.0-m_om)))/(1.0-m_om);
	case OL_2:
		return m_t_h*(1.0-m_om/sqrt(m_om-1.0)*atan(sqrt(m_om-1.0)))/(1.0-m_om);
	case OL_3:
		return 2.0*m_t_h/3.0;
		break;
	case OM:
		return m_t_h*asinh(1.0/(sqrt(1.0/m_ol-1.0)))/sqrt(m_ol);
		break;
	case A1:
		return ta1(0);
		break;
	case A2_1:
	case A2_2:
		return ta2(0);
		break;
	case OM_OL_1:
		return 2.0*m_t_h*asinh(sqrt(m_ol/m_om))/(3.0*sqrt(m_ol));
		break;
	}
	return -1;
}

double flrw::age(double z) const {
	switch (m_case) {
	case OM_OL_0:
		return m_t_h/(1+z);
		break;
	case OL_1:
	case OL_2:
	case OL_3:
		return tolz(z);
		break;
	case OM:
		return tomz(z);
		break;
	case A1:
		return ta1(z);
		break;
	case A2_1:
	case A2_2:
		return ta2(z);
		break;
	case OM_OL_1:
		return tb(z);
		break;
	}
	return -1;
}

// m_ol=0 CASE: OL_1,OL_2,OL_3
double flrw::tolz(double z) const {
	const double pre0=1.0-m_om;
	const double prez=sqrt(1.0+m_om*z);
	switch (m_case) {
	case OL_1:
		//ol=0 0<m_om<1
		return m_t_h*(prez/(1.0+z)-m_om/sqrt(pre0)*atanh(sqrt(pre0)/prez))/pre0;
	case OL_2:
		//m_ol=0 m_om>1
		return m_t_h*(prez/(1.0+z)-m_om/sqrt(-pre0)*atan(sqrt(-pre0)/prez))
				/pre0;
	case OL_3:
		// m_ol=0 m_om=1
		return 2.0*m_t_h/(3.0*(1.0+z)*sqrt(1.0+z));
	}
	return -1.;
}

// m_om=0 CASE: OM
double flrw::tomz(double z) const {
	return m_t_h*asinh(1.0/((1.0+z)*sqrt(1.0/m_ol-1.0)))/sqrt(m_ol);
}

// CASE A1
double flrw::ta1(double z) const {
	gsl_set_error_handler_off();
	int status=0;
	gsl_sf_result result;
	const double vk=pow(m_kap*(m_b-1.0)+sqrt(m_b*(m_b-2.0)), 1.0/3.0);
	const double y1=(m_kap*(vk+1.0/vk)-1.0)/3.0;
	const double A=sqrt(y1*(3.0*y1+2.0));
	const double k=sqrt((2.0*A+m_kap*(1+3.0*y1))/(4.0*A));
	double arg0=m_kap*y1+m_om*(1.0+z)/fabs(m_ok);
	double phi=acos((arg0-A)/(arg0+A));
	const double sin_phi=sin(phi);
	double n=-0.25*(A+m_kap*y1)*(A+m_kap*y1)/(A*m_kap*y1);
	double phi_null=asin(1/sqrt(-n));
	double arg2, arg3;
	double hm, hp;
	double pre;
	double arg1=(1.0+z)*m_om/m_ok;
	if (fabs(1.0+n*sin_phi*sin_phi)<EPS) {
		n=-y1*(1.0+y1)/(A-m_kap*y1)/(A-m_kap*y1);
		//  EQUATION 22, a very special case of b = 27*(2+sqrt(2))/8.
		if (fabs(1.0+n*sin_phi*sin_phi)<EPS) {
			phi=acos(-1.0-arg1/M_SQRT2+1.0-arg1);
			const double arg2=(1.0-arg1)*sqrt(arg1*arg1+(arg1+M_SQRT1_2)*(1
					+M_SQRT1_2));
			const double arg3=(M_SQRT2-1.0)*(arg1+M_SQRT2+1.0)*sqrt((1.0
					+M_SQRT1_2)*(M_SQRT1_2-arg1));
			return 0.25*m_t_h/sqrt(m_ol)*((M_SQRT2-1.0)*gsl_sf_ellint_F(phi,
					0.5*sqrt(1+2*M_SQRT2), PREC)+log(fabs((arg2+arg3)/(arg2
					-arg3))));
		}
		//        EQUATION 8.
		else {
			// En mi paper, 7
			arg2=arg1*arg1*(1.0+arg1);
			arg2=sqrt((1.0+y1)*((1.0+y1)*y1*y1-arg2));
			arg2*=2*m_kap*y1;
			arg1*=arg1*(A-m_kap*y1)-2*m_kap*(1.0+y1)*y1*y1;
			hm=arg1-arg2;
			hp=arg1+arg2;
			arg1=gsl_sf_ellint_F(phi, k, PREC)/(m_kap*y1*sqrt(A));
			status=gsl_sf_ellint_P_e(phi, n, k, PREC, &result);
			if (status) {
				//cerr<<"error: "<<gsl_strerror(status)<<" on "
				//    <<__FILE__<<" "<<__LINE__<<endl;
				return integrate_time(z);
			}
			arg2=(A-m_kap)/(y1*(1.0+y1)*sqrt(A))*result.val;
			arg3=log(fabs(hm/hp))/(m_kap*y1*sqrt(m_kap*(y1+1)));
			pre=0.5*m_om/fabs(m_ok)/sqrt(fabs(m_ok));
			return m_t_h*pre*(arg1+arg2+arg3);
		}
	}
	//       EQUATION 10.
	else {
		// En mi paper, 9
		hm=sqrt(((1.0+y1)*(y1-arg1))/(y1*y1+(1.0+arg1)*(y1+arg1)));
		arg1=-gsl_sf_ellint_F(phi, k, PREC)/(A+m_kap*y1);
		status=gsl_sf_ellint_P_e(phi, n, k, PREC, &result);
		if (status) {
			//cerr<<"error: "<<gsl_strerror(status)<<" on "
			//    <<__FILE__<<" "<<__LINE__<<endl;
			return integrate_time(z);
		}
		arg2=-0.5*(A-m_kap*y1)/(m_kap*y1*(A+m_kap*y1))*result.val;
		arg3=-0.5*(sqrt(A/(m_kap*(y1+1.0)))/(m_kap*y1))*log(fabs((1.0-hm)/(1.0
				+hm)));
		return m_t_h*m_om/(m_ok*sqrt(A*fabs(m_ok)))*(arg1+arg2+arg3);
	}
	return -1.0;
}

double flrw::ta2(double z) const {
	//       EQUATION 19, a very special case of b = 2.
	if (m_case==A2_1) {
		const double arg=(1.0+z)*m_om/m_ok;
		return -m_t_h*(M_SQRT3*log(1.0-2.0/(sqrt(1.0/3.0-arg)+1.0))+log(1.0+2.0/(sqrt(1.0-3.0*arg)-1.0)))/sqrt(m_ol);
	}
	//       EQUATION 15.
	if (m_case==A2_2) {
		int status=0;
		gsl_sf_result result;
		const double arg=acos(1.0-m_b)/3.0;
		const double arg1=cos(arg)/3.0;
		const double arg2=sin(arg)/M_SQRT3;
		const double y1=-1.0/3.0+arg1+arg2;
		const double y2=-1.0/3.0-2.0*arg1;
		const double arg3=y1-y2;
		const double y3=y1-2.0*arg2;
		const double phi=asin(sqrt(arg3/(y1-m_om*(1.0+z)/m_ok)));
		const double k=sqrt((y1-y3)/arg3);
		const double n=-y1/arg3;
		const double arg4=y1*fabs(m_ok)*sqrt(-arg3*m_ok);
		status=gsl_sf_ellint_P_e(phi, n, k, PREC, &result);
		if (status) {
			//cerr<<"error: "<<gsl_strerror(status)<<" on "
			//    <<__FILE__<<" "<<__LINE__<<endl;
			return integrate_time(z);
		}
		return 2.0*m_t_h*m_om/arg4*(result.val-gsl_sf_ellint_F(phi, k, PREC));
	}
	return -1.0;
}

double flrw::tb(double z) const {
	// m_om+m_ol=1
	return 2.0*m_t_h*asinh(sqrt(m_ol/(m_om*(1.0+z*(3.0+z*(3.0+z))))))/(3.0
			*sqrt(m_ol));
}

double flrw::integrate_time(double z) const {
	gsl_integration_workspace* w=gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_function F;
	F.function=&helper_fun_time;
	F.params=static_cast<void*>(const_cast<flrw*>(this));
	const int status=gsl_integration_qagiu(&F, z, 0, 1e-7, 1000, w, &result,
			&error);
	if (status) {
		std::cerr<<"error: "<<gsl_strerror(status)<<" on " <<__FILE__<<" "<<__LINE__<<std::endl;
    }
    return m_t_h*result;
  }
} //namespace metrics
} // namespace milia
