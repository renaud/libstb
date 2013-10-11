#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <time.h>

gsl_rng *gslr = NULL;
int verbose = 0;

#define DIM 8
#define TERMS 4
#define BCYCLES 50000

double pvec[DIM];
double qvec[DIM];

void yap(char *s) {
  fprintf(stderr,"quit: %s\n", s);
  exit(1);
}

int main(int argc, char* argv[])
{
  int i, c, d, iter, ITER=0;
  int seed=0;
  double alpha = 0.0;
  double beta = 0.5;
  double q1=0, q1sq=0, q2=0, q2sq=0;

  /*
   *  default values for args
   */
  ITER = 5000;
  while ( (c=getopt(argc, argv,"a:b:vs:I:"))>=0 ) {
    switch ( c ) {
    case 'v':     
      verbose++;
      break;    
    case 'b':
      if ( !optarg || sscanf(optarg,"%lf",&beta)!=1 )
        yap("Need a valid 'b' argument\n");
      break;
    case 'a':
      if ( !optarg || sscanf(optarg,"%lf",&alpha)!=1 )
        yap("Need a valid 'a' argument\n");
      break;
    case 'I':
      if ( !optarg || sscanf(optarg,"%d",&ITER)!=1 )
        yap("Need a valid 'I' argument\n");
      break;
    case 's':
      if ( !optarg || sscanf(optarg,"%d",&seed)!=1 )
        yap("Need a valid 's' argument\n");
      break;
    }
  }

  /*
   *   set random number generator
   */
  if ( seed ) {
    gsl_rng_default_seed = seed;
  } else {
    gsl_rng_default_seed = time(NULL);
  }
  gslr = gsl_rng_alloc(gsl_rng_rand48);
  if ( !gslr ) 
    yap("cannot allocate gsl_rng\n");
  fprintf(stderr, "Setting seed = %lu\n", gsl_rng_default_seed);

  /*
   *   make p
   */
  { 
    double tot = 0;
    for (i=0; i<DIM; i++) {
      tot += pvec[i] = DIM-i;
    }
    for (i=0; i<DIM; i++) {
      pvec[i] /= tot;
    }
  }
  
  for (iter=1; iter<=ITER; iter++) {
    double tot;
    double prod;
    /*
     *  generate q via Dirichlet
     */
    tot = 0;
    for (d=0; d<DIM; d++) {
      tot += qvec[d] = gsl_ran_gamma(gslr, beta*pvec[d], 1);
    }
    for (d=0; d<DIM; d++) {
      qvec[d] /= tot;
    }
    if ( verbose ) {
      printf("q1 = ");
      for (d=0; d<DIM; d++)
	printf(" %lg", qvec[d]);
      printf("\n");
    }
    prod = 1.0;
    for (d=0; d<TERMS; d++)
      prod *= qvec[d];    
    q1 += prod;
    q1sq += prod*prod;
    /*
     *  generate q via DP
     */
    for (d=0; d<DIM; d++)
      qvec[d] = 0;
    tot = 1.0;
    for (i=1; i<=BCYCLES; i++) {
      double v = gsl_ran_beta(gslr, 1-alpha, beta+i*alpha);
      double p = tot*v;
      unsigned int n;
      double flat = gsl_ran_flat(gslr,0,1);
      tot *= (1-v);
      for (n=0; n<DIM && flat>0; n++)
	flat -= pvec[n];
      assert(n>0);
      qvec[n-1] += p;
    }
    if ( verbose ) {
      printf("q2 = ");
      for (d=0; d<DIM; d++)
	printf(" %lg", qvec[d]);
      printf("\n");
    }
    prod = 1.0;
    for (d=0; d<TERMS; d++)
      prod *= qvec[d];    
    q2 += prod;
    q2sq += prod*prod;
    if ( (iter%100)==0 && iter>0 ) {
      printf("mean(%d) = %lg +/- %lg, %lg+/- %lg\n", iter,
	     q1/iter, sqrt(q1sq/iter-(q1/iter)*(q1/iter))/sqrt(iter-1),
	     q2/iter, sqrt(q2sq/iter-(q2/iter)*(q2/iter))/sqrt(iter-1) );
    }
  }
  printf("mean(%d) = %lg +/- %lg, %lg+/- %lg\n", iter,
	 q1/iter, sqrt(q1sq/iter-(q1/iter)*(q1/iter))/sqrt(iter-1),
	 q2/iter, sqrt(q2sq/iter-(q2/iter)*(q2/iter))/sqrt(iter-1) );
  return 1;
}
