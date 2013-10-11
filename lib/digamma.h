/*
 *  Define digamma() and polygamma() here.
 *
 *  We use the GSL special functions, GSL is available
 *  on most platforms including MAC and Windows.
 *
 *  Another option is the Cephes library with psi() and zeta()
 *           http://www.netlib.org/cephes/
 *  which have nice, small self contained functions without all the
 *  GSL cruft.
 */

#ifndef __DIGAMMA_H
#define __DIGAMMA_H

#include <gsl/gsl_sf.h>

#define digamma(x) gsl_sf_psi(x)
#define polygamma(i,x) gsl_sf_psi_n(i,x)

/*
 *  this is given in digamma.c
 */
double digammaInv(double x);

#endif

