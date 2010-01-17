/*
 * Copyright 2009-2010 Sergio Pascual
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

#include "FlrwAge.h"
#include "milia/flrw_nat.h"
#include "milia/exception.h"

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_math.h>

namespace
{
  const double ITOL = 1.0e-7;
  const size_t NPOINTS = 1000;
  struct Params
  {
      double m;
      double v;
  };

  double helper_fun_time(double z, void* pars)
  {
    Params* pmetric = static_cast<Params*> (pars);
    const double om = pmetric->m;
    const double ol = pmetric->v;
    return 1. / ((1. + z) * (sqrt(gsl_pow_2 (1. + z) * (1. + om * z) - z * ol
        * (2. + z))));
  }

  double integrate_time(double z, double m, double v)
  {
    gsl_error_handler_t* oldhandler = gsl_set_error_handler_off();
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(NPOINTS);
    double result, error;
    gsl_function F;
    F.function = &helper_fun_time;
    Params p = { m, v };
    F.params = static_cast<void*> (&p);
    const int status = gsl_integration_qagiu(&F, z, 0, ITOL, NPOINTS, w,
        &result, &error);
    gsl_integration_workspace_free(w);
    gsl_set_error_handler(oldhandler);
    if (status)
    {
      throw milia::exception(std::string("gsl error: ") + gsl_strerror(status));
    }
    return result;
  }
}
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(FlrwAge);

void FlrwAge::setUp()
{
}

void FlrwAge::tearDown()
{
}

void FlrwAge::testAge()
{
  const double z = 1;
  for (double om = 0.00; om < 2.01; om += 0.01)
    for (double ol = 0.00; ol < 2.01; ol += 0.01)
    {
      bool infll = false;
      double iv = -1;
      double cv = -1;
      try
      {
        iv = integrate_time(z, om, ol);
      } catch (milia::exception& e)
      {
        infll = true;
      }
      try
      {
        const milia::metrics::flrw_nat m(om, ol);
        cv = m.age(z);
      } catch (milia::recollapse& e)
      {
        infll = true;
      } catch (milia::no_big_bang& e)
      {
        infll = true;
      }
      if (not infll)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(iv, cv, ITOL);
      }
    }
}
