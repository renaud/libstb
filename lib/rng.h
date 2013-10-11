/*
 *    a simple interface to the chosen RNG
 */
#ifndef __RNG_H
#define __RNG_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef gsl_rng *rngp_t;

#define rng_unit(rng) gsl_ran_flat(rng, 0, 1)
#define rng_beta(rng, a, b) gsl_ran_beta(rng, a, b)
#define rng_gamma(rng, a) gsl_ran_gamma(rng, a, 1)
#define rng_gaussian(rng, a) gsl_ran_gaussian(rng, a)
#endif
