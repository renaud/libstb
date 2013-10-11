/*
 * Simple slice sampler on log concave distribution
 * Copyright (C) 2012 Wray Buntine 
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *  Simple slice sampler on log-concave distribution
 *          
 */
#include <assert.h>
#include <math.h>

#include "rng.h"
// #include "yaps.h"

#define TOOMANY 200
/*
 *    on input, *xp should be an approximate maxima
 *               to make sampling behave beter
 *    post() is a log posterior calc, so handle accordingly
 *    do a few loops to let burn in
 *
 *   returns 1 on error
 */
int SliceSimple(double *xp, 
		double (*post)(double, void *), 
		double *bounds,
		rngp_t rng, int loops, void *pars) {
  double x = *xp, y;
  double range[2];
  int tries;
  if ( x<bounds[0] || x>bounds[1] ) {
    fprintf(stderr,
	    "SliceSimple: input value %lf outside bounds [%lg,%lg]\n",
	    x, bounds[0], bounds[1]);
    return 1;
  }
  //  yaps_message("SliceSimple: x=%lg\n", x);
  while ( loops-->0 ) {
    y = post(x,pars);
    range[0] = bounds[0];
    range[1] = bounds[1];
    // yaps_message(" start p(%lg)=%lg, range=[%lg,%lg]\n",x,y,range[0],range[1]);
    y += log(rng_unit(rng));
    // yaps_message("  y=%lg\n", y);
    for (tries=1; tries<TOOMANY; tries++) {
      x = range[0] + rng_unit(rng)*(range[1]-range[0]);
      if ( post(x,pars)>y ) {
	*xp = x;
	// yaps_message(" got p(%lg)=%lg after %d tries, range=[%lg,%lg]\n",
	//	    x, post(x,pars), tries, range[0], range[1]);
	break;
      }
      /*
       *  reduce bounds, assumes posterior is unimodal
       */
      if ( x<*xp ) 
	range[0] = x;
      else
	range[1] = x;
      // yaps_message("   p(%lg)=%lg, r=[%lg,%lg]\n",
      //             x,post(x,pars),range[0],range[1]);
    }
    if (tries>=TOOMANY ) {
      fprintf(stderr,
	      "SliceSimple: giving up after %d tries, range=[%lg,%lg]\n",
	      TOOMANY, range[0], range[1]);
      return 1;
    }
  }
  return 0;
}
