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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/cbrt.hpp>

#include "flrw.h"

using std::sqrt;
using std::atan;
using std::log;
using std::abs;

using boost::math::ellint_1;
using boost::math::ellint_3;
using boost::math::asinh;
using boost::math::atanh;
using boost::math::cbrt;

namespace {
const double M_SQRT3 = sqrt(3);
const double M_4THRT3 = sqrt(M_SQRT3);
}

namespace milia {
namespace metrics {

double flrw::age() const {
	switch (m_case) {
	case OM_OV_0:
		return m_t_h;
		break;
	case OV_1:
		return m_t_h * (1 - m_om / sqrt(1 - m_om) * atanh(sqrt(1 - m_om))) / (1
				- m_om);
	case OV_2:
		return m_t_h * (1 - m_om / sqrt(m_om - 1) * atan(sqrt(m_om - 1))) / (1
				- m_om);
	case OV_3:
		return 2 * m_t_h / 3;
		break;
	case OM:
		return m_t_h * asinh(1 / (sqrt(1 / m_ov - 1))) / sqrt(m_ov);
		break;
	case A1:
		return ta1(0);
		break;
	case A2_1:
	case A2_2:
		return ta2(0);
		break;
	case OM_OV_1:
		return 2 * m_t_h * asinh(sqrt(m_ov / m_om)) / (3 * sqrt(m_ov));
		break;
	}
	return -1;
}

double flrw::age(double z) const {
	switch (m_case) {
	case OM_OV_0:
		return m_t_h / (1 + z);
		break;
	case OV_1:
	case OV_2:
	case OV_3:
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
	case OM_OV_1:
		return tb(z);
		break;
	}
	return -1;
}

// ol=0 CASE: OL_1,OL_2,OL_3
double flrw::tolz(double z) const {
	const double pre0 = 1 - m_om;
	const double prez = sqrt(1 + m_om * z);
	switch (m_case) {
	case OV_1:
		//ol=0 0<om<1
		return m_t_h * (prez / (1 + z) - m_om / sqrt(pre0) * atanh(sqrt(pre0)
				/ prez)) / pre0;
	case OV_2:
		//ol=0 om>1
		return m_t_h * (prez / (1 + z) - m_om / sqrt(-pre0) * atan(sqrt(-pre0)
				/ prez)) / pre0;
	case OV_3:
		// ol=0 om=1
		return 2 * m_t_h / (3 * (1 + z) * sqrt(1 + z));
	}
	return -1.;
}

// om=0 CASE: OM
double flrw::tomz(double z) const {
	return m_t_h * asinh(1 / ((1 + z) * sqrt(1 / m_ov - 1))) / sqrt(m_ov);
}

// CASE A1
double flrw::ta1(double z) const {
	const double vk = cbrt(m_kap * (m_b - 1) + sqrt(m_b * (m_b - 2)));
	const double y1 = (m_kap * (vk + 1 / vk) - 1) / 3.;
	const double A = sqrt(y1 * (3 * y1 + 2));
	const double k = sqrt((2 * A + m_kap * (1 + 3 * y1)) / (4 * A));
	double arg0 = m_kap * y1 + m_om * (1 + z) / abs(m_ok);
	double phi = acos((arg0 - A) / (arg0 + A));
	const double sin_phi = sin(phi);
	double n = -0.25 * (A + m_kap * y1) * (A + m_kap * y1) / (A * m_kap * y1);
	double arg2, arg3;
	double hm, hp;
	double pre;
	double arg1 = (1 + z) * m_om / m_ok;
	if (1 + n * sin_phi * sin_phi == 0) {
		n = -y1 * (1 + y1) / (A - m_kap * y1) / (A - m_kap * y1);
		//  EQUATION 22, a very special case of b = 27*(2+sqrt(2))/8.
		if (1 + n * sin_phi * sin_phi == 0) {
			phi = acos(-1 - arg1 / M_SQRT2 + 1 - arg1);
			const double arg2 = (1 - arg1) * sqrt(arg1 * arg1 + (arg1
					+ M_SQRT1_2) * (1 + M_SQRT1_2));
			const double arg3 = (M_SQRT2 - 1) * (arg1 + M_SQRT2 + 1) * sqrt((1
					+ M_SQRT1_2) * (M_SQRT1_2 - arg1));
			return 0.25 * m_t_h / sqrt(m_ov) * ((M_SQRT2 - 1) * ellint_1(0.5
					* sqrt(1 + 2 * M_SQRT2), phi) + log(abs((arg2 + arg3)
					/ (arg2 - arg3))));
		}
		//        EQUATION 8.
		else {
			arg2 = arg1 * arg1 * (1 + arg1);
			arg2 = sqrt((1 + y1) * ((1 + y1) * y1 * y1 - arg2));
			arg2 *= 2 * m_kap * y1;
			arg1 *= arg1 * (A - m_kap * y1) - 2 * m_kap * (1 + y1) * y1 * y1;
			hm = arg1 - arg2;
			hp = arg1 + arg2;
			arg1 = ellint_1(k, phi) / (m_kap * y1 * sqrt(A));
			arg2 = (A - m_kap) / (y1 * (1 + y1) * sqrt(A))
					* ellint_3(k, n, phi);
			arg3 = log(abs(hm / hp)) / (m_kap * y1 * sqrt(m_kap * (y1 + 1)));
			pre = 0.5 * m_om / abs(m_ok) / sqrt(abs(m_ok));
			return m_t_h * pre * (arg1 + arg2 + arg3);
		}
	}
	//       EQUATION 10.
	else {
		hm = sqrt(((1 + y1) * (y1 - arg1)) / (y1 * y1 + (1 + arg1)
				* (y1 + arg1)));
		arg1 = -ellint_1(k, phi) / (A + m_kap * y1);
		arg2 = -0.5 * (A - m_kap * y1) / (m_kap * y1 * (A + m_kap * y1))
				* ellint_3(k, n, phi);
		arg3 = -0.5 * (sqrt(A / (m_kap * (y1 + 1))) / (m_kap * y1)) * log(abs(
				(1.0 - hm) / (1.0 + hm)));
		return m_t_h * m_om / (m_ok * sqrt(A * abs(m_ok))) * (arg1 + arg2
				+ arg3);
	}
	return -1.0;
}

double flrw::ta2(double z) const {
	//       EQUATION 19, a very special case of b = 2.
	if (m_case == A2_1) {
		const double arg = (1 + z) * m_om / m_ok;
		return -m_t_h * (M_SQRT3 * log(1 - 2 / (sqrt(1. / 3. - arg) + 1))
				+ log(1 + 2 / (sqrt(1 - 3 * arg) - 1))) / sqrt(m_ov);
	}
	//       EQUATION 15.
	if (m_case == A2_2) {
		const double arg = acos(1 - m_b) / 3;
		const double arg1 = cos(arg) / 3;
		const double arg2 = sin(arg) / M_SQRT3;
		const double y1 = -1. / 3. + arg1 + arg2;
		const double y2 = -1. / 3. - 2 * arg1;
		const double arg3 = y1 - y2;
		const double y3 = y1 - 2 * arg2;
		const double phi = asin(sqrt(arg3 / (y1 - m_om * (1 + z) / m_ok)));
		const double k = sqrt((y1 - y3) / arg3);
		const double n = -y1 / arg3;
		const double arg4 = y1 * abs(m_ok) * sqrt(-arg3 * m_ok);
		return 2 * m_t_h * m_om / arg4 * (ellint_3(k, n, phi)
				- ellint_1(k, phi));
	}
	return -1.0;
}

double flrw::tb(double z) const {
	// om+ol=1
	return 2 * m_t_h * asinh(sqrt(m_ov / (m_om * (1 + z * (3 + z * (3 + z))))))
			/ (3 * sqrt(m_ov));
}
} // namespace metrics
} // namespace milia
