/*
 * digamma routines
 * Copyright (C) 2010 Wray Buntine 
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
 */

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "digamma.h"

/*
 *    Minka's algorithm
 *    gets overflow on values > 80;
 *    large negative values return near zero;
 */
double digammaInv(double x) {
  double guess;
  int i;
  if ( x < -2.22 )
    guess = -1 / (x - digamma(1.0));
  else
    guess = exp(x) + 0.5;
  for (i=0; i<5; i++ ) {
    guess -= (digamma(guess) - x)/ polygamma(1,guess);
  }
  return guess;
}
