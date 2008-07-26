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

// $Id $

#ifndef _MILIA_METRIC_H_
#define _MILIA_METRIC_H_

namespace milia {

/**
 * The Friedmann-Lemaître-Robertson-Walker metric
 */
class metric {
public:
	metric(double hubble_constant, double matter_density, double lambda_density);
private:
};

/**
 * Long name of the metric
 */
typedef metric flrw_metric;

} // namespace milia


#endif /* _MILIA_METRIC_H_ */