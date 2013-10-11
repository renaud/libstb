/*
 * Stirling Number table handling, extra features
 * Copyright (C) 2009-2010 Wray Buntine 
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
 *     
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "yaps.h"
#include "sxtable.h"
#include "sapprox.h"

#ifndef NDEBUG
/*
 *   define TESTS to create a second much larger SNM table
 *   against which the thinned one and its computations are checked;
 *   see SX_safe() ;
 *   report of errors done at end
 */
// #define TESTS
/*
 *   does more extensive and safer version of recSNM()
 */
// #define TEST_RECSNM
#endif

extern int verbose;

#ifdef TESTS
static int abserr(double v1, double v2) {
  if ( fabs(v1-v2)/v1 > 0.0001 )
    return 1;
  return 0;
}
#endif

/*
 *   stores the S table as a float
 */
static float **SX_m;
static float **dSX_m;
static float SX_a;
/*
 *  but all calculations done as a double, so we need to
 *  keep the N and M frontier in double
 */
static double *SX_Nfrontier;
static double *SX_Mfrontier;
static double *dSX_Nfrontier;
static double *dSX_Mfrontier;
static double *SX_Mfrontier_old;
/*
 *  dimensions
 */
static int usedM, usedN, fullM;
static unsigned int memory, start_memory;
/*
 *   if non-zero, we only keep SNM(*,M) for every skip-th M,
 *   after the initial M<=fullM
 */
static int init = 1;
static int skip = 1;
static int fix;
static int diva;

/*
 *   buffer of size 2 to save a copy
 */
static struct SX_store_s {
   float **SX_m;
   float **dSX_m;
   float SX_a;
   double *SX_Nfrontier, *SX_Mfrontier, *SX_Mfrontier_old;
   double *dSX_Nfrontier, *dSX_Mfrontier;
   int usedM, usedN, fullM;
   unsigned int memory, start_memory;
   int skip, fix, diva;
   int primary;
 } SX_store[2];

// #define REPORT_SNM
#ifdef REPORT_SNM
static FILE *stfp = NULL;
#endif

static double logadd(double V, double lp) {
  if ( lp>V ) {
    // swap so V is bigger
    double t = lp;
    lp = V;
    V = t;
  }
  return V + log(1.0+exp(lp-V));
}

#ifdef TESTS
static float **SX_m_test;
static int usedM_test, usedN_test;
static long test_total = 0;
static long test_wrong = 0;
static long test_nocheck = 0;

static void makeSX_main(int maxN, int maxM, float a, int skipin,
		       int fixin, int diva);
void SX_make(int maxN, int maxM, float a, int skipin, int fixin, int diva) {
  /*
   *    call once to create a larger table for checking answers
   */
  makeSX_main(maxN*3, maxM*4, a, skipin, fixin, diva);
  usedM_test = usedM;
  usedN_test = usedN;
  yaps_message("SX_make test build SNM for N<%d and M<=%d\n",
	     usedN, usedM);
  SX_m_test = SX_m;
  /*
   *    now this is the version used
   */
  SX_m = NULL;
  usedM = usedN = 0;
  makeSX_main(maxN, maxM, a, skipin, fixin, diva);
}

/*
 *  fills S table
 */
static void makeSX_main(int maxN, int maxM, float a, int skipin, int fixin, int diva) {
#else
/*
 *  fills S table
 */
  void SX_make(int maxN, int maxM, float a, int skipin, int fixin, int divin) {
#endif
  int N, M;

  if ( verbose )
    yaps_message("Making S for N=%d/%d M=%d a=%f with div=%d\n", 
	       maxN, skipin, maxM, a, divin);

  if ( init ) {
    /*
     *   tag these as unused
     */
    SX_store[0].primary = -1;
    SX_store[1].primary = -1;
    init = 0;
  }
  usedN = maxN;
  skip = skipin;
  fix = fixin;
  diva = divin;
  SX_a = a;
  if ( skip<=1 )
    skip = 1;
  else {
    int diff = maxN%skip;
    if ( diff>0 )
      maxN += skip-diff;
    assert(maxN%skip==0);
  }
  fullM = usedM = maxM;
  
  SX_m = malloc(sizeof(SX_m[0])*(maxM+1));
  memory = maxM*maxN;
  if ( !SX_m )
    yaps_quit("Cannot allocate SX_m\n");
  SX_m[0] = NULL;
  for (M=1; M<=maxM; M++) {
    SX_m[M] = malloc(sizeof(SX_m[0][0])*(maxN+1));
    if ( !SX_m[M] )
      yaps_quit("Cannot allocate SX_m[%d]\n", M);
  }

  SX_Nfrontier = malloc(sizeof(SX_Nfrontier[0])*maxN);
  if ( !SX_Nfrontier )
    yaps_quit("Cannot allocate SX_Nfrontier\n");
  SX_Mfrontier = malloc(sizeof(SX_Mfrontier[0])*(maxM+1));
  if ( !SX_Mfrontier )
    yaps_quit("Cannot allocate SX_Mfrontier\n");
  memory += (maxN+maxM+1)*2;
  start_memory = memory;
  /*
   *  all values outside bounds to 0
   */
  for (N=1; N<maxN; N++) {
    for (M=N+1; M<=maxM; M++)
      SX_m[M][N] =  -HUGE_VAL;
  }

  if ( diva ) {
    dSX_m = malloc(sizeof(SX_m[0])*(maxM+1));
    memory = maxM*maxN;
    if ( !dSX_m )
      yaps_quit("Cannot allocate dSX_m\n");
    dSX_m[0] = NULL;
    for (M=1; M<=maxM; M++) {
      dSX_m[M] = malloc(sizeof(SX_m[0][0])*maxN);
      if ( !dSX_m[M] )
	yaps_quit("Cannot allocate dSX_m[%d]\n", M);
    }
    
    dSX_Nfrontier = malloc(sizeof(SX_Nfrontier[0])*maxN);
    if ( !dSX_Nfrontier )
      yaps_quit("Cannot allocate dSX_Nfrontier\n");
    dSX_Mfrontier = malloc(sizeof(SX_Mfrontier[0])*(maxM+1));
    SX_Mfrontier_old = malloc(sizeof(SX_Mfrontier[0])*(maxM+1));
    if ( !SX_Mfrontier_old )
      yaps_quit("Cannot allocate dSX_Mfrontier\n");
    memory += (maxN+2*(maxM+1))*2;
    start_memory = memory;
    /*
     *  all values outside bounds to 0
     */
    for (N=1; N<maxN; N++) {
      for (M=N+1; M<=maxM; M++)
	dSX_m[M][N] =  -HUGE_VAL;
    }
  }
  
#ifdef REPORT_SNM
  stfp = fopen("st.txt","w");
#endif
  
  SX_m[1][1] = 0;
  /*
   *   on entry, SX_Mfrontier[M] holds the value for SNM(1,M)
   */
  SX_Mfrontier[1] = 0;
  if ( diva ) {
    dSX_m[1][1] = 0;
    /*
     *   on entry, dSX_Mfrontier[M] holds the value for dSNM(1,M)
     */
    dSX_Mfrontier[1] = 0;
  }
  for (N=2; N<maxN; N++) {
    /*
     *  we overwrite frontier (which has SNM(N-1,*) stored),
     *  need to do in reverse order so insitu is OK
     */
    for (M=maxM; M>1; M--) {
      if ( diva ) 
	SX_Mfrontier_old[M] = SX_Mfrontier[M];
      if ( M>N )
	continue;
      if ( N==M )
	SX_Mfrontier[M] = 0;
      else
	SX_Mfrontier[M] = logadd(SX_Mfrontier[M-1], 
				log(N-(M*a)-1.0) + SX_Mfrontier[M]);
      // SX_m[M][N] = logadd(SX_m[M-1][N-1],log(N-(M*a)-1.0)+SX_m[M][N-1]);
      if ( !isfinite(SX_Mfrontier[M]) ) {
	fprintf(stderr," SX_NM(%d,%d) gone inf. during adding\n", N, M);
	exit(1);
      } 
    }
    if ( diva ) 
      SX_Mfrontier_old[1] = SX_Mfrontier[1];
    SX_Mfrontier[1] = log(N-a-1.0) + SX_Mfrontier[1];
    for (M=1; M<=maxM && M<=N; M++) {
      SX_m[M][N] = SX_Mfrontier[M];
    }
    if ( N>=maxM )
      SX_Nfrontier[N] = SX_Mfrontier[maxM];
    
    if ( diva ) {
      /*
       *  we overwrite frontier (which has dSNM(N-1,*) stored),
       *  need to do in reverse order so insitu is OK
       */
      for (M=maxM; M>1; M--) {
	if ( M>N )
	  continue;
	if ( N==M )
	  dSX_Mfrontier[M] = 0;
	else
	  dSX_Mfrontier[M] = dSX_Mfrontier[M-1] * exp(SX_Mfrontier_old[M-1]-SX_Mfrontier[M])
	    + ((N-(M*a)-1.0)*dSX_Mfrontier[M] - M) * exp(SX_Mfrontier_old[M]-SX_Mfrontier[M]);
	if ( !isfinite(dSX_Mfrontier[M]) ) {
	  fprintf(stderr," dSX_NM(%d,%d) gone inf. during adding\n", N, M);
	  exit(1);
	} 
      }
      dSX_Mfrontier[1] = -1/(N-a-1.0) + dSX_Mfrontier[1];
      for (M=1; M<=maxM && M<=N; M++) {
	dSX_m[M][N] = dSX_Mfrontier[M];
      }
      if ( N>=maxM )
	dSX_Nfrontier[N] = dSX_Mfrontier[maxM];
    }
  }
  
#ifdef CHECK_SNM
  /*
   *   on exit, SX_Mfrontier[0:maxM] holds the value for SNM(usedN-1,0:maxM)
   *      SX_Nfrontier[maxM:maxN-1] holds the value for SNM(maxM:usedN-1,maxM)
   */
  { int bad=0;
   for (N=2; N<maxN; N++) {
     for (M=1; M<=4; M++) {
       double v1 = SX_m[M][N];
       double v2 = S_approx(N,M,a);
       if ( fabs(v1-v2)/v2>0.0001 ) {
	 yaps_message("Bad SNM(%d,%d) %lf vs approx=%lf\n",
		     N,M,v1,v2);
	 bad++;
       }
     }
     if ( bad>200 ) 
       yaps_quit("give up!\n");
   }
  }
#endif
}

static void extendT(int maxM) {
  int N, M;
  int lastM = usedM;

  if ( fix )
    yaps_quit("Asking to extend T, but is set to fix\n");

  if ( verbose )
	yaps_message("Extending M for S from %d to %d\n", usedM, maxM);
  usedM = maxM;

  SX_m = realloc(SX_m, sizeof(SX_m[0])*(maxM+1));
  if ( !SX_m )
    yaps_quit("Cannot allocate SX_m\n");
  for (M=lastM+1; M<=maxM; M++) {
    SX_m[M] = malloc(sizeof(SX_m[0][0])*(usedN/skip+1));
    if ( !SX_m[M] )
      yaps_quit("Cannot allocate SX_m[%d]\n", M);
  }
  memory += (usedN/skip+1)*(maxM-lastM);

  SX_Mfrontier = realloc(SX_Mfrontier, sizeof(SX_Mfrontier[0])*(maxM+1));
  if ( !SX_Mfrontier )
    yaps_quit("Cannot allocate SX_Mfrontier\n");
  memory += 2*(maxM-lastM);

  /*
   *  all values outside bounds to 0
   */
  for (N=1; N<maxM; N++) {
    if ( skip>1 && N%skip != 0 )
      continue;
    M = N+1;
    if ( M<=lastM )
      M = lastM+1;
    for ( ; M<=maxM; M++) {
      SX_m[M][N] =  -HUGE_VAL;
    }
  }
  for (M=lastM+1; M<=maxM; M++) {
    double this_lastM = SX_Nfrontier[M];
    SX_Nfrontier[M] = 0;
    for (N=M+1; N<usedN; N++) {
      double save_lastM = SX_Nfrontier[N];   //  save SNM(N,M-1)
      SX_Nfrontier[N] = logadd(this_lastM, log(N-(M*SX_a)-1.0) + SX_Nfrontier[N-1]);
      // SX_m[M][N] = logadd(SX_m[M-1][N-1], log(N-(M*a)-1.0) + SX_m[M][N-1]);
      this_lastM = save_lastM;
      if ( !isfinite(SX_Nfrontier[N]) ) {
	fprintf(stderr," SX_NM(%d,%d) gone inf. during adding\n", N, M);
	exit(1);
      }
    }
    for (N=M; N<usedN; N++) {
      if ( skip>1 ) {
	if ( N%skip != 0 )
	  continue;
	SX_m[M][N/skip] = SX_Nfrontier[N];
      } else
	SX_m[M][N] = SX_Nfrontier[N];
    }
    if ( M<usedN )
      SX_Mfrontier[M] = SX_Nfrontier[usedN-1];
  }
  // yaps_message("Extended M for S\n");
}

static void extendN(int maxN) {
  int N, M;
  int lastN = usedN;

  if ( verbose )
    yaps_message("Extending N for S from %d to %d\n", usedN, maxN);
  usedN = maxN;

  for (M=1; M<=fullM; M++) {
    SX_m[M] = realloc(SX_m[M], sizeof(SX_m[0][0])*maxN);
    if ( !SX_m[M] )
      yaps_quit("Cannot allocate SX_Nfrontier\n");
  }
  memory += (maxN-lastN)*(fullM);
  for (M=fullM+1; M<=usedM; M++) {
    SX_m[M] = realloc(SX_m[M], sizeof(SX_m[0][0])*(maxN/skip+1));
    if ( !SX_m[M] )
      yaps_quit("Cannot allocate SX_Nfrontier\n");
  }
  memory += (usedM-fullM)*(maxN/skip-lastN/skip);

  SX_Nfrontier = realloc(SX_Nfrontier, sizeof(SX_Nfrontier[0])*maxN);
  if ( !SX_Nfrontier )
    yaps_quit("Cannot allocate SX_Nfrontier\n");
  memory += (maxN-lastN)*2;

  for (N=lastN; N<usedN; N++) {
    for (M=usedM; M>1; M--) {
      if ( M>N )
 	continue;
      if ( N==M )
	SX_Mfrontier[M] = 0;
      else
	SX_Mfrontier[M] = logadd(SX_Mfrontier[M-1], log(N-(M*SX_a)-1.0) + SX_Mfrontier[M]);
      // SX_m[M][N] = logadd(SX_m[M-1][N-1], log(N-(M*SX_a)-1.0) + SX_m[M][N-1]);
      if ( !isfinite(SX_Mfrontier[M]) ) {
	fprintf(stderr," SX_NM(%d,%d) gone inf. during adding\n", N, M);
	exit(1);
      }
    }
    SX_Mfrontier[1] = log(N-SX_a-1.0) + SX_Mfrontier[1];
    for (M=1; M<=usedM && M<=N; M++) {
      if ( skip>1 && M>fullM ) {
	if ( M%skip!=0 ) 
	  continue;
	SX_m[M][N/skip] = SX_Mfrontier[M];
      } else {
	SX_m[M][N] = SX_Mfrontier[M];
      }
    }
    if ( N>=usedM )
      SX_Nfrontier[N] = SX_Mfrontier[usedM];
  }
  // yaps_message("Extended N for S\n");
}

void SX_check(int N, int T) {
  if ( N+1 >= usedN ) 
    extendN(N*1.1);
  if ( !fix &&  T+1 >= usedM ) {
    int ext = T*1.1;
    if ( T<200 )
      ext = T+20;
    extendT(ext);
  }
}

static double recSNM_da(int N, int T) {
  if ( T<=fullM )
    return dSX_m[T][N];
  if ( N%skip==0 ) {
    return dSX_m[T][N/skip];
  } 
  if ( N==T )
    return 0;
  yaps_quit("recSNM_da() not done\n");
  return 0;
}
static double recSNM(int N, int T) {
  if ( T<=fullM )
    return SX_m[T][N];
  if ( N%skip==0 ) {
    return SX_m[T][N/skip];
  } 
  if ( N==T )
    return 0;
  return logadd(recSNM(N-1,T-1), log(N-(T*SX_a)-1.0) + recSNM(N-1,T));
}
#ifdef TEST_RECSNM
static double recSNM_print(int N, int T, float a) {
  if ( T<=fullM ) {
    yaps_message("recSNM(%d,%d<=%d) = %f",
		N, T, fullM, SX_m[T][N]);
#ifdef TESTS
    if ( abserr(SX_m[T][N], SX_m_test[T][N]) )
      yaps_message(", error!");
#endif
    yaps_message("\n");
    return SX_m[T][N];
  } 
  if ( N%skip==0 ) {
    double ts =  SX_m[T][N/skip];
    yaps_message("recSNM(%d/%d,%d) = %f",
		N, skip, T, ts );
#ifdef TESTS
    if ( abserr(ts, SX_m_test[T][N]) )
      yaps_message(", error!");
#endif
    yaps_message("\n");
    return ts;
  } 
  if ( N==T )
    return 0;
  {
    double rec1 = recSNM(N-1,T-1,SX_a);
    double rec2 = recSNM(N-1,T,SX_a);
    assert(N-(T*SX_a)-1.0>0);
    yaps_message("recSNM(%d,%d) = la(%lf, %lf+%lf) = %lf",
		N, T, rec1, log(N-(T*SX_a)-1.0), rec2,
		logadd(rec1,log(N-(T*SX_a)-1.0) + rec2) );
#ifdef TESTS
    if ( abserr(rec1, SX_m_test[T-1][N-1]) )
      yaps_message(", first value error!");
    if ( abserr(rec2, SX_m_test[T][N-1]) )
      yaps_message(", second value error!");
#endif
    yaps_message("\n");
    return logadd(rec1,log(N-(T*SX_a)-1.0) + rec2);
  }
}
#endif

double SX_safe_da(int N, int T) {
  if ( N+1 >= usedN && T<=4 ) 
    return S_approx_da(N, T, SX_a);
  SX_check(N, T);
  if ( skip>1 && (T>fullM) )
    return recSNM_da(N,T);
  if ( N==T )
    return 0;
  // assert( T<usedM && N<usedN );
  return dSX_m[T][N];
}

double SX_safe(int N, int T) {
  if ( N+1 >= usedN && T<=4 ) 
    return S_approx(N, T, SX_a);
  SX_check(N, T);
#ifdef REPORT_SNM
  if ( T>20 ) {
    fprintf(stfp, "%d %d\n", N, T);
    fflush(stfp);
  }
#endif
#ifdef TESTS
  {
    double val;
    if ( skip>1 && T>fullM )
      val = recSNM(N,T);
    else
      val = SX_m[T][N];
    if ( N<usedN_test && T<usedM_test ) {
      test_total++;
      if ( (val-SX_m_test[T][N])/val > 0.00001 ) {
	yaps_message("SNM(%d,%d) computed=%lf precomputed=%lf\n",
		    N, T, val, SX_m_test[T][N]);
	recSNM_print(N,T);
	test_wrong++;
	if ( test_wrong>10 )
	  yaps_quit("Quiting\n");
      }
    } else
      test_nocheck++;
    return val;
  }
#else
  if ( skip>1 && (T>fullM) )
    return recSNM(N,T);
  if ( N==T )
    return 0;
  return SX_m[T][N];
#endif
}

void SX_report(void) {
  yaps_message("SNM initial memory = %u, extended = %u floats\n", 
	     start_memory, memory);
#ifdef TESTS
  if ( test_nocheck ) 
    yaps_message(", but %ld untested", test_nocheck);
#endif
  yaps_message("\n");
}


#define saves(I,name) SX_store[I].name = name
#define restores(I,name) name = SX_store[I].name

static void freeSS(int I) {
  int M;
  //yaps_message("freeSS(%d): ", I);
  if ( SX_store[I].primary<0 ) {
    //yaps_message("not allocated\n");
    return;
  }
  //yaps_message("usedN=%d, usedM=%d, fullM=%d, primary=%d\n",
  //	      SX_store[I].usedN, SX_store[I].usedM, 
  //	      SX_store[I].fullM, SX_store[I].primary);
  free(SX_store[I].SX_Nfrontier);
  free(SX_store[I].SX_Mfrontier);
  if ( SX_store[I].diva ) {
    free(SX_store[I].dSX_Nfrontier);
    free(SX_store[I].dSX_Mfrontier);
    free(SX_store[I].SX_Mfrontier_old);
  }
  skip = 1;
  for (M=1; M<=SX_store[I].usedM; M++)
    free(SX_store[I].SX_m[M]);
  free(SX_store[I].SX_m);
  if ( SX_store[I].diva ) {
    for (M=1; M<=SX_store[I].usedM; M++)
      free(SX_store[I].dSX_m[M]);
    free(SX_store[I].dSX_m);
    SX_store[I].dSX_m = NULL;
  }
  SX_store[I].SX_m = NULL;
  SX_store[I].SX_Nfrontier = NULL;
  SX_store[I].SX_Mfrontier = NULL;
  SX_store[I].primary = -1;
}

void SX_save(int I, int primary) {
  freeSS(I);
  SX_store[I].primary = primary;
   saves(I,SX_m);
   saves(I,SX_a);
   saves(I,dSX_m);
   saves(I,SX_Nfrontier);
   saves(I,SX_Mfrontier);
   saves(I,SX_Mfrontier_old);
   saves(I,dSX_Nfrontier);
   saves(I,dSX_Mfrontier);
   saves(I,usedM);
   saves(I,usedN);
   saves(I,fullM);
   saves(I,memory);
   saves(I,start_memory);
   saves(I,skip);
   saves(I,fix);
   saves(I,diva);
 }
 
static void restoreS(int I) {
  assert(SX_store[I].primary>=0);
   restores(I,SX_m);
   restores(I,SX_a);
   restores(I,dSX_m);
   restores(I,SX_Nfrontier);
   restores(I,SX_Mfrontier);
   restores(I,SX_Mfrontier_old);
   restores(I,dSX_Nfrontier);
   restores(I,dSX_Mfrontier);
   restores(I,usedM);
   restores(I,usedN);
   restores(I,fullM);
   restores(I,memory);
   restores(I,start_memory);
   restores(I,skip);
   restores(I,fix);
   restores(I,diva);
  }

 /*
  *     use parameters in the store to rebuild,
  *     but leave store alone
  */
 void SX_rebuild(float apar, int primary) {
   if ( primary ) {
     /*
      *   get and store a new primary table, so
      *   must have diva set;
      *   discard any non-primary table if space needed
      */
     if ( SX_store[0].primary>=0
	  && SX_store[0].diva
	  && apar==SX_store[0].SX_a ) {
       restoreS(0);
       SX_store[0].primary = 1;
       if ( SX_store[1].primary==1 )
	 SX_store[1].primary = 0;
     } else if ( SX_store[1].primary>=0
		 && SX_store[1].diva
		 && apar==SX_store[1].SX_a ) {
       restoreS(1);
       SX_store[1].primary = 1;
       if ( SX_store[0].primary==1 )
	 SX_store[0].primary = 0;
     } else if ( SX_store[0].primary<1 ) {
       if ( SX_store[1].primary==1 ) 
	 SX_store[1].primary = 0;
       SX_make(usedN, fullM, apar, skip, fix, 1);
       if ( fullM<usedM )
	 extendT(usedM);
       SX_save(0,1);
     } else {
       assert(SX_store[1].primary<1 );
       if ( SX_store[0].primary==1 ) 
	 SX_store[0].primary = 0;
       SX_make(usedN, fullM, apar, skip, fix, 1);
       if ( fullM<usedM )
	 extendT(usedM);
       SX_save(1,1);
     }
   } else {
     /*
      *   get and store a new non-primary table
      *   replace any non-primary table if space needed
      */
     if ( SX_store[0].primary>=0 
	  && apar==SX_store[0].SX_a ) {
       restoreS(0);
     } else if ( SX_store[1].primary>=0 
		 && apar==SX_store[1].SX_a ) {
       restoreS(1);
     } else if ( SX_store[0].primary<1 ) {
       SX_make(usedN, fullM, apar, skip, fix, 0);  //WLB should be 0
       if ( fullM<usedM )
	 extendT(usedM);
       SX_save(0,0);
     } else {
       assert(SX_store[1].primary<1 );
       SX_make(usedN, fullM, apar, skip, fix, 0);  //WLB should be 0
       if ( fullM<usedM )
	 extendT(usedM);
       SX_save(1,0);
     }
   }
 }
   

void SX_free(void) {
#ifdef TESTS
  yaps_message("SNM testing:  %ld wrong out of %ld", test_wrong, test_total);
#endif
  freeSS(0);
  freeSS(1);
}

float SX_alpha(void) {
  return SX_a;
}
