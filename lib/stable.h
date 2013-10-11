/*
 * Stirling Number table handling for Pitman-Yor processing
 * Copyright (C) 2009-2012 Wray Buntine 
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
 *   This is a modified version of earlier code
 *   with a rebuilt standardised interface.
 *   For definitions and math. details, see
 *   	http://arxiv.org/abs/1007.0296
 *
 */

#ifndef __STABLE_H
#define __STABLE_H

#include "stdint.h"

/*
 *   values for flags, once set are fixed
 *   S_STABLE - build and maintain the S table
 *   S_UVTABLE - build and maintain the U+V+S1 tables
 *   S_FLOAT - keep final table values in float, but all 
 *             intermediate calcs done in double
 *   S_VERBOSE - print occasional stats to stdout
 *   S_QUITONBOUND - if the 2 maximum bounds overran, then die,
 *                   otherwise return 0 or log(0)
 */
#define S_STABLE 1
#define S_UVTABLE 2
#define S_FLOAT 4
#define S_VERBOSE 8
#define S_QUITONBOUND 16

/*
 *  dimensions, parameters, data stored here;
 *  should be "private", so don't look at this;
 *  some pointers are set to NULL if unused
 */
typedef struct stable_s {
  /*
   *    inclusive bounds fixed for life
   */
  unsigned maxM, maxN;    
  /*
   *      inclusive current bounds, may be increased automatically
   *      as needed
   */
  unsigned usedM, usedN;
  /*
   *      used for deallocation, initial (usedM-2),
   *      because S[0:startMM][.] is single block of memory
   */
  unsigned startM;        
  /*
   *   stores  S[n][m] = log S^{n+3}_{m+2,a}  for n,m >=0
   */
  double **S;
  float **Sf;
  double *SfrontN;    //  frontier, S[n][m] for n=usedN-3, m>=2
  double *SfrontM;    //  frontier, S[n][m] for m=usedM-2, n>usedM-3
  /*
   *   stores  V[n][m] = V^{n+2}_{m+2,a}  for n,m >=0
   */
  double **V;
  float **Vf;
  double *VfrontN;   //  frontier, V[n][m] for n=usedN-2, m>=2
  double *VfrontM;   //  frontier, V[n][m] for m=usedM-2, n>usedM-2
  /*
   *   stores  S1[n] = log S^{n+1}_{1,a}, for  0 <= n < usedN
   */  
  unsigned usedN1;        // has memory allocated for usedN1>=usedN
  double *S1;
   /*
   *   stores  log gamma(1-a)
   */   
  double lga;
  double a;
  /*
   *    flags
   */
  uint32_t flags;
  uint32_t memalloced;
  char *tag;   /*  name */
} stable_t;


/*
 *  fills S table with initN and initM values;
 *     - these are reset to be 10 if any less,
 *     - and be maxN/M if anymore
 *  flags and maxN/M bounds are fixed for life;
 *  
 *  S_FLOAT ::  all larger tables are stored as floats, though
 *              intermediate calculations are all double,
 *              halves memory requirements
 *
 *  return NULL on error
 */
stable_t *S_make(unsigned initN, unsigned initM, unsigned maxN, unsigned maxM, 
		 double a, uint32_t flags);
void S_tag(stable_t *S, char *tag);

/*
 *  fill with new values for different "a"
 *  return non-zero on error
 */
int S_remake(stable_t *sp, double a);

void S_free(stable_t *sp);

/*
 *  all the table value routines instigate a rebuild
 *  if current tables not big enough, and required values
 *  are within maxN and maxM bounds, otherwise return 0;
 *  caller needs to handle the "over bounds" case themselves
 */

/*
 *   return  log S^n_{m,a}
 *   i.e.,  does bounds checking and extends table if needed
 *   return log(0) = -INFINITY if over maxN/M bounds, or no table
 */
double S_S(stable_t *sp, unsigned n, unsigned m);

/*
 *   return log S^n_{1,a} = log (gamma(n-a)/gamma(1-a) )
 *   for n >= 1
 *   using cache if possible
 *   return log(0) = -INFINITY if over maxN bounds, or no table
 */
double S_S1(stable_t *sp, unsigned n);

/*
 *    return U^n_{m,a} = S^{n+1}_{m,a}/S^{n}_{m,a} 
 *    computed from T_V(), for m>=1
 *    return 0 if illegal inputs or out of bounds, or no table
 */
double S_U(stable_t *sp, unsigned n, unsigned m);
/*
 *    return S_U(...)*S_V(...) using one less call
 */
double S_UV(stable_t *sp, unsigned n, unsigned m);
    
/*
 *   return V^n_{m,a} = S^{n}_{m,a}/S^{n}_{m-1,a} 
 *   from tables, for m>=2
 *   return 0 if illegal inputs or out of bounds
 */
double S_V(stable_t *sp, unsigned n, unsigned m);

/*
 *    print a simple one line of stats to *fp
 */
void S_report(stable_t *sp, FILE *fp);

#endif
