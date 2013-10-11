/*
 * Stirling Number table handling
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
 *   This is a rewrite of various earlier implementations with a different
 *   data layout and a few more options.
 *
 *   Note the code has a float and a double version.  The float version
 *   keeps two double vectors representating the current boundary of
 *   the Stirling number matrix.  This way it can be extended in either
 *   direction without sacrificing precision.  However, it makes
 *   extensions somewhat complicated, because one needs to keep
 *   track of the boundary during computation.
 *     
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "stable.h"
#include "yaps.h"


#define MEMALLOCED
#ifdef MEMALLOCED
/********************************************************************
 *    this is a fudge to allow recording of memory allocation sizes;
 *    we store the allocation size along with the memory
 */
//   guess at basic memory requirements of malloc()
#define memsize(s) ((((s)+9)/8)*8)

void *malloc_hook(stable_t *sp, size_t size) {
  void *ptr;
  size_t ms = memsize(size);
  size += sizeof (size_t);
  ptr = malloc(size);
  if ( !ptr ) return NULL;
  *(size_t *) ptr = size;
  if ( sp ) sp->memalloced += ms;
  return ((size_t *) ptr) + 1;
}

void *realloc_hook(stable_t *sp, void *ptr, size_t size) {
  size_t ms = memsize(size);
  size += sizeof (size_t);
  ptr = realloc((void *) (((size_t *) ptr) - 1),size);
  if ( !ptr ) return NULL;
  *(size_t *) ptr = size;
  if ( sp ) sp->memalloced += ms;
   return ((size_t *) ptr) + 1;
}

void free_hook (stable_t *sp, void *ptr) {
  if ( sp ) sp->memalloced -= memsize(* (((size_t *) ptr) - 1));
  ptr = (void *) (((size_t *) ptr) - 1);
  free(ptr);
}

#define realloc(p,x) realloc_hook(sp,p,x)
#define malloc(x) malloc_hook(sp,x)
#define free(x) free_hook(sp,x)
/********************************************************************/
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

void S_tag(stable_t *S, char *tag) {
  S->tag = tag;
}

stable_t *S_make(unsigned initN, unsigned initM, unsigned maxN, unsigned maxM, 
		 double a, uint32_t flags) {
  int N;
  stable_t *sp = NULL;
  sp = malloc(sizeof(stable_t));
  if ( !sp ) 
    return NULL;
  
  if ( maxM<10 )
    maxM = 10;
  if ( maxN<maxM )
    maxN = maxM;
  if ( initM<10 )
    initM = 10;
  if ( initN<initM )
    initN = initM;
  if ( initN>maxN )
    initN = maxM;
  if ( initN>maxN )
    initN = maxN;

  if ( (flags&S_STABLE)==0 && (flags&S_UVTABLE)==0 )
    return NULL;
  
  sp->tag = NULL;
  sp->memalloced = 0;
  sp->flags = flags;
  sp->maxN = maxN;
  sp->maxM = maxM;
  sp->usedN = initN;
  sp->usedM = initM;
  sp->usedN1 = initN;
  sp->startM = initM;
  sp->S = NULL;
  sp->SfrontN = sp->SfrontM = sp->S1 = NULL;
  sp->V = NULL;
  sp->VfrontN = sp->VfrontM = NULL;
  sp->Sf = sp->Vf = NULL;
  
  sp->S1 = malloc(sizeof(sp->S1[0])*(initN));
  if ( !sp->S1 ) {
    free(sp);
    return NULL;
  }
  if ( flags&S_STABLE ) {
    if ( flags&S_FLOAT ) {
      /*
       *  allocate frontier
       */
      sp->SfrontN = malloc(sizeof(sp->SfrontN[0])*(initM-1));
      if ( !sp->SfrontN ) {
	S_free(sp);  	return NULL;
      }
      /*
       *    sets diagonal entry of S since the loop writing
       *    SfrontN never does the diagnal itself
       */
      memset(sp->SfrontN,0,sizeof(sp->SfrontN[0])*(initM-1));
      sp->SfrontM = malloc(sizeof(sp->SfrontM[0])*(initN-initM+1));
      if ( !sp->SfrontM ) {
	S_free(sp);  	return NULL;
      }
      /*
       *  allocate sp->Sf[] as vector of vectors
       */
      sp->Sf = malloc(sizeof(sp->Sf[0])*(sp->usedN-2));
      if ( !sp->Sf ) {
	S_free(sp);  	return NULL;
      }
      memset(sp->Sf,0,sizeof(sp->Sf[0])*(sp->usedN-2));
      /*
       *   allocate sp->Sf[0][.] to sp->Sf[startM-3][.] in one block
       */
      sp->Sf[0] = malloc(sizeof(sp->Sf[0][0])*(sp->startM-1)*(sp->startM-2)/2);
      if ( !sp->Sf[0] ) {
	S_free(sp);  return NULL;
      }
      for (N=1; N<=sp->startM-3; N++) 
	sp->Sf[N] = sp->Sf[N-1] + N;
      /*
       *   allocate remaining sp->Sf[N][.] as vectors
       */
      assert(sp->startM-2+sp->Sf[sp->startM-3]-sp->Sf[0]==
	     (sp->startM-1)*(sp->startM-2)/2);
      for (N=sp->startM-2; N<=sp->usedN-3; N++) {
	sp->Sf[N] = malloc(sizeof(sp->Sf[0][0])*(sp->usedM-1));
	if ( !sp->Sf[N] ) {
	  S_free(sp);  return NULL;
	}
      }  
    } else { 
      sp->S = malloc(sizeof(sp->S[0])*sp->usedN);
      if ( !sp->S ) {
	S_free(sp);  return NULL;
      }
      /*
       *   allocate sp->S[0][.] to sp->S[startM-3][.] in one block
       */
      sp->S[0] = malloc(sizeof(sp->S[0][0])*(sp->startM-1)*(sp->startM-2)/2);
      if ( !sp->S[0] ) {
	S_free(sp);  return NULL;
      }
      for (N=1; N<=sp->startM-3; N++) 
	sp->S[N] = sp->S[N-1] + N;
      /*
       *   allocate remaining sp->S[N][.] as vectors for N>=usedM+1
       *   which store values M=2,...,usedM, so need (usedM-1) space
       */
      assert(sp->startM-2+sp->S[sp->startM-3]-sp->S[0]==
	     (sp->startM-1)*(sp->startM-2)/2);
      for (N=sp->startM-2; N<=sp->usedN-3; N++) {
	sp->S[N] = malloc(sizeof(sp->S[0][0])*(sp->usedM-1));
	if ( !sp->S[N] ) {
	  S_free(sp);  return NULL;
	}
      }  
    }
  }
  if ( flags&S_UVTABLE ) {
    if ( flags&S_FLOAT ) {
     /*
       *  allocate frontier
       */
      sp->VfrontN = malloc(sizeof(sp->VfrontN[0])*(initM-1));
      if ( !sp->VfrontN ) {
	S_free(sp);  	return NULL;
      }
      /*
       *    sets diagonal entry of V since the loop writing
       *    VfrontN never does the diagnal itself
       */
      memset(sp->VfrontN,0,sizeof(sp->VfrontN[0])*(initM-1));
      sp->VfrontM = malloc(sizeof(sp->VfrontM[0])*(initN-initM+1));
      if ( !sp->VfrontM ) {
	S_free(sp);  	return NULL;
      }
      
      sp->Vf = malloc(sizeof(sp->Vf[0])*sp->usedN);
      if ( !sp->Vf ) {
	S_free(sp);  return NULL;
      }
      /*
       *   allocate sp->Vf[0][.] to sp->Vf[startM-2][.] in one block
       */
      sp->Vf[0] = malloc(sizeof(sp->Vf[0][0])*(sp->startM-1)*(sp->startM)/2);
      if ( !sp->Vf[0] ) {
	S_free(sp);  return NULL;
      }
      for (N=1; N<=sp->startM-2; N++) 
	sp->Vf[N] = sp->Vf[N-1] + N;
      /*
       *   allocate remaining sp->Vf[N][.] as vectors for N>=usedM+1
       *   which store values M=2,...,usedM, so need (usedM-1) space
       */
      assert(sp->startM-1+sp->Vf[sp->startM-2]-sp->Vf[0]==
	     (sp->startM-1)*(sp->startM)/2);
      for (N=sp->startM-1; N<=sp->usedN-2; N++) {
	sp->Vf[N] = malloc(sizeof(sp->Vf[0][0])*(sp->usedM-1));
	if ( !sp->Vf[N] ) {
	  S_free(sp);  return NULL;
	}
      }  
    } else {
      sp->V = malloc(sizeof(sp->V[0])*sp->usedN);
      if ( !sp->V ) {
	S_free(sp);  return NULL;
      }
      /*
       *   allocate sp->V[0][.] to sp->V[startM-2][.] in one block
       */
      sp->V[0] = malloc(sizeof(sp->V[0][0])*(sp->startM-1)*(sp->startM)/2);
      if ( !sp->V[0] ) {
	S_free(sp);  return NULL;
      }
      for (N=1; N<=sp->startM-2; N++) 
	sp->V[N] = sp->V[N-1] + N;
      /*
       *   allocate remaining sp->V[N][.] as vectors for N>=usedM+1
       *   which store values M=2,...,usedM, so need (usedM-1) space
       */
      assert(sp->startM-1+sp->V[sp->startM-2]-sp->V[0]==
	     (sp->startM-1)*(sp->startM)/2);
      for (N=sp->startM-1; N<=sp->usedN-2; N++) {
	sp->V[N] = malloc(sizeof(sp->V[0][0])*(sp->usedM-1));
	if ( !sp->V[N] ) {
	  S_free(sp);  return NULL;
	}
      }  
    }
  }

  /*
   *  this is where we actually build the Stirling numbers
   */
  S_remake(sp,a);
  return sp;
}


/*
 *    assumes s->usedN/M already set to new values and memory filled
 *    startN/M = 0  --->  refill everything
 *    startN/M > 0  --->  memory extended so refill from here,
 *                        i.e.,  these were *last* values set, start +1
 */
static int S_remake_part(stable_t *sp, double a, 
			 unsigned startN, unsigned startM) {
  int N, M;
  
  if ( startN==0 )
    startM = 0;
  sp->a = a;
  sp->lga = lgamma(1.0-a);

  // yaps_message("S_remake_part(a=%lf,N=%u, M=%u)\n", a, startN, startM);

  /*
   *  need to reset at sp->S1[] least to sp->usedN1;
   *  up to sp->usedN used by sp->S[][],
   *  and data needs overwriting up to sp->usedN1
   */
  if ( startN==0 ) {
    sp->S1[0] = 0;
    N = 2;
  } else {
    N = startN+1;
    assert(sp->S1[startN-1]>0 );
  }
  for ( ; N<=sp->usedN; N++)
    sp->S1[N-1] = sp->S1[N-2] + log(N-1-a);

  if ( startN==0 )
    //   a has changed, so reset others
    for ( ; N<=sp->usedN1; N++)
      sp->S1[N-1] = 0;
  
  if ( sp->flags&S_STABLE ) {
    if ( (sp->flags&S_FLOAT)==0 ) {
      if ( startM>0 && startM<sp->usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+1; N<=startN; N++) {
	  for (M=startM+1; M<N && M<=sp->usedM; M++) {
	    sp->S[N-3][M-2] = 
	      logadd(log(N-M*a-1.0)+((M<N-1)?sp->S[N-4][M-2]:0), 
				     sp->S[N-4][M-3]);
	    assert(isfinite(sp->S[N-3][M-2]));
	  }
	}
      }
      /*
       *   now fill from 0 after startN ... just like usual
       */
      if ( startN==0 ) {
	sp->S[0][0] = logadd(sp->S1[1],log(2-2*a));
	N = 4;
      } else {
	N = startN+1;
	assert(sp->S[N-4][0]>0);
      }
      for (; N<=sp->usedN; N++) {
	sp->S[N-3][0] = logadd(log(N-2*a-1.0)+sp->S[N-4][0], 
			       sp->S1[N-2]);
	for (M=3; M<=sp->usedM && M<N; M++) {
	  sp->S[N-3][M-2] = 
	    logadd(log(N-M*a-1.0)+((M<N-1)?sp->S[N-4][M-2]:0), sp->S[N-4][M-3]);
	  assert(isfinite(sp->S[N-3][M-2]));
	}
      }
    } else {
      /*
       *   computation done in double by storing in sp->SfrontN+M[]
       */
      if ( startM>0 && startM<sp->usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+2; N<=startN; N++) {
	  double lastS;
	  if (startM+1<N-1)
	    lastS = sp->Sf[N-4][startM-1];
	  else
	    lastS = 0 ;
	  sp->Sf[N-3][startM-1] = sp->SfrontN[startM-1]
	    = logadd(log(N-(startM+1)*a-1.0)+lastS,
		     sp->SfrontM[N-startM-2]);
	  for (M=startM+2; M<N && M<=sp->usedM; M++) {
	    double saveS = sp->SfrontN[M-2];
	    if ( M==N-1 ) saveS = 0;
	    sp->SfrontN[M-2] = logadd(log(N-M*a-1.0)+saveS, lastS);
	    sp->Sf[N-3][M-2] = sp->SfrontN[M-2];
	    assert(isfinite(sp->Sf[N-3][M-2]));
	    lastS = saveS;
	  }
	  //   save the SfrontM value
	  if ( N>sp->usedM )
	    sp->SfrontM[N-sp->usedM-1] = sp->SfrontN[sp->usedM-2];
	}
      }
      if ( startN==0 ) {
	sp->Sf[0][0] = sp->SfrontN[0] = logadd(sp->S1[1],log(2-2*a));
	N = 4;
      } else {
	N = startN+1;
	assert(sp->Sf[N-4][0]>0);
      }
      for ( ; N<=sp->usedN; N++) {
	double lastS;
	lastS = sp->SfrontN[0];
	sp->Sf[N-3][0] = sp->SfrontN[0] =
	  logadd(log(N-2*a-1.0)+lastS, sp->S1[N-2]);
	for (M=3; M<=sp->usedM && M<N; M++) {
	  double saveS = sp->SfrontN[M-2];
	  if (M==N-1) saveS = 0;
	  sp->SfrontN[M-2] =
	    logadd(log(N-M*a-1.0)+saveS, lastS);
	  sp->Sf[N-3][M-2] = sp->SfrontN[M-2];
	  assert(isfinite(sp->Sf[N-3][M-2]));
	  lastS = saveS;
	}
	//   save the SfrontM value
	if ( N>sp->usedM )
	  sp->SfrontM[N-sp->usedM-1] = sp->SfrontN[sp->usedM-2];
      }
    }
  }    
  if ( sp->flags&S_UVTABLE ) {
    if ( (sp->flags&S_FLOAT)==0 ) {
      if ( startM>0 && startM<sp->usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+1; N<=startN; N++) {
	  for (M=startM+1; M<=N && M<=sp->usedM; M++) {
	    sp->V[N-2][M-2] = 
	      (1.0+((M<N)?((N-1-M*a)*sp->V[N-3][M-2]):0))
	      / (1.0/sp->V[N-3][M-3]+(N-1-(M-1)*a));
	  }
	}
      }
      /*
       *   now fill from 0 after startN ... just like usual
       */
      if ( startN==0 ) {
	sp->V[0][0] = 1.0/(1.0-a);
	N = 3;
      } else {
	N = startN+1;
	assert(sp->V[N-3][0]>0);
      }
      for (; N<=sp->usedN; N++) {
	sp->V[N-2][0] = (1.0+(N-1-2*a)*sp->V[N-3][0])/(N-1-a);
	for (M=3; M<=sp->usedM && M<=N; M++) {
	  sp->V[N-2][M-2] = 
	    (1.0+((M<N)?((N-1-M*a)*sp->V[N-3][M-2]):0))
	    / (1.0/sp->V[N-3][M-3]+(N-1-(M-1)*a));
	}
      }
    } else {
      if ( startM>0 && startM<sp->usedM ) {
	/*
	 *   extend for M upto startN
	 */
	for (N=startM+1; N<=startN; N++) {
	  double lastS;
	  if ( startM+1<N) 
	    lastS = sp->VfrontN[startM-1];
	  else 
	    lastS = 0;  
	  sp->Vf[N-2][startM-1] = sp->VfrontN[startM-1] = 
	    (1.0+(N-1-(startM+1)*a)*lastS)
	    / (1.0/sp->VfrontM[N-1-startM]+(N-1-(startM)*a));
	  for (M=startM+2; M<=N && M<=sp->usedM; M++) {
	    double saveS = sp->VfrontN[M-2];
	    assert(lastS!=0);
	    sp->Vf[N-2][M-2] = sp->VfrontN[M-2] = 
	      (1.0+((M<N)?((N-1-M*a)*saveS):0))
	      / (1.0/lastS+(N-1-(M-1)*a));
	    lastS = saveS;
	  }
	  //   save the VfrontM value
	  if ( N>=sp->usedM )
	    sp->VfrontM[N-sp->usedM] = sp->VfrontN[sp->usedM-2];
 	}
      }
      /*
       *   now fill from 0 after startN ... just like usual
       */
      if ( startN==0 ) {
	sp->Vf[0][0] = sp->VfrontN[0] = 1.0/(1.0-a);
	N = 3;
      } else {
	N = startN+1;
	assert(sp->Vf[N-3][0]>0);
      }
      for (; N<=sp->usedN; N++) {
	double lastS;
	lastS = sp->VfrontN[0];
	sp->Vf[N-2][0] = sp->VfrontN[0] =
	  (1.0+(N-1-2*a)*lastS)/(N-1-a);
	for (M=3; M<=sp->usedM && M<=N; M++) {
	  double saveS = sp->VfrontN[M-2];
	  assert(lastS!=0);
	  sp->Vf[N-2][M-2] = sp->VfrontN[M-2] = 
	    (1.0+((M<N)?((N-1-M*a)*saveS):0))
	    / (1.0/lastS+(N-1-(M-1)*a));
	  lastS = saveS;
	}
	//   save the VfrontM value
	if ( N>=sp->usedM )
	  sp->VfrontM[N-sp->usedM] = sp->VfrontN[sp->usedM-2];
      }
    }
  }
  return 0;
}

int S_remake(stable_t *sp, double a) {
  int ret = S_remake_part(sp, a, 0, 0);
  if ( !ret && (sp->flags&S_VERBOSE) ) 
    S_report(sp, stderr);
  return ret;
}

/*
 *    assume tables filled to usedN,usedM;
 *    request (N,M) table value;
 *    fiddle maxM, maxN to make it non-trivial
 *    create extra space first;
 *    then remake table values;
 *    return non-zero on error
 */
static int S_extend(stable_t *sp, int N, int M) {
  int n;
  /*
   *  N shouldn't be too big or small
   */
  if ( N<sp->usedN )
    N = sp->usedN;
  if ( N>sp->maxN )
    N = sp->maxN;

  /*
   *   N increase should not be trivial
   */
  if ( N>sp->usedN ) {
    if ( N<sp->usedN*1.1 )
      N = sp->usedN*1.1;
    if ( N<sp->usedN+50 )
      N = sp->usedN+50;
    /*
     *  reset if made too big
     */
    if ( N>sp->maxN ) {
      N = sp->maxN;
    }
  }

  /*
   *  M shouldn't be too big or small
   */
  if ( M<sp->usedM )
    M = sp->usedM;
  if ( N<M )
    M = N;
  if ( M>sp->maxM )
    M = sp->maxM;

  /*
   *   M increase should not be trivial
   */
  if ( M>sp->usedM ) {
    if ( M<sp->usedM*1.1 )
      M = sp->usedM*1.1;
    if ( M<sp->usedM+50 )
      M = sp->usedM+50;
    /*
     *  reset if made too big
     */
    if ( M>sp->maxM ) {
      M = sp->maxM;
    }
    if ( M>sp->usedN ) {
      M = sp->usedN;
    }
  }
  /*
   *   N and M values now set
   */

  if ( M>sp->usedM ) {
    /*
     *    extend size of existing vectors;
     *    note S/Sf/V/Vf are triangular up to n=usedM, so ignore that
     */
    for (n=sp->usedM+1; n<=sp->usedN; n++) {
      if ( sp->flags&S_STABLE ) {
	if ( (sp->flags&S_FLOAT)==0 ) {
	  sp->S[n-3] = realloc(sp->S[n-3],sizeof(sp->S[0][0])*(M-1));
	  if ( !sp->S[n-3] ) {
	    S_free(sp);
	    return 1;
	  }
	} else {
	  sp->Sf[n-3] = realloc(sp->Sf[n-3],sizeof(sp->Sf[0][0])*(M-1));
	  if ( !sp->Sf[n-3] ) {
	    S_free(sp);
	    return 1;
	  }
	}
      }
      if ( sp->flags&S_UVTABLE ) {
	if ( (sp->flags&S_FLOAT)==0 ) {
	  sp->V[n-2] = realloc(sp->V[n-2],sizeof(sp->V[0][0])*(M-1));
	  if ( !sp->V[n-2] ) {
	    S_free(sp);
	    return 1;
	  }
	} else {
	  sp->Vf[n-2] = realloc(sp->Vf[n-2],sizeof(sp->Vf[0][0])*(M-1));
	  if ( !sp->Vf[n-2] ) {
	    S_free(sp);
	    return 1;
	  }
	}
      }
    }
    if ( sp->flags&S_STABLE && sp->flags&S_FLOAT ) {
      sp->SfrontN = realloc(sp->SfrontN,sizeof(sp->SfrontN[0])*(M-1));
      if ( !sp->SfrontN ) {
	S_free(sp);
	return 1;
      }
    }
    if ( sp->flags&S_UVTABLE && sp->flags&S_FLOAT ) {
      sp->VfrontN = realloc(sp->VfrontN,sizeof(sp->VfrontN[0])*(M-1));
      if ( !sp->VfrontN ) {
	S_free(sp);
	return 1;
      }
    }
  }
  /*
   *    now create new vectors in S/Sf/V/Vf
   */
  if ( N>sp->usedN ) {
    /*
     *   extend size of S/Sf/SfrontN
     */
    if ( sp->flags&S_STABLE ) {
      sp->S1 = realloc(sp->S1,sizeof(sp->S1[0])*N);
      if ( !sp->S1 ) {
	S_free(sp);
	return 1;
      }  
      if ( (sp->flags&S_FLOAT)==0 ) {
	sp->S = realloc(sp->S,sizeof(sp->S[0])*N);
	if ( !sp->S ) {
	  S_free(sp);
	  return 1;
	}
      } else {
	sp->Sf = realloc(sp->Sf,sizeof(sp->Sf[0])*N);
	if ( !sp->Sf ) {
	  S_free(sp);
	  return 1;
	}
	sp->SfrontM = realloc(sp->SfrontM,sizeof(sp->SfrontM[0])*(N-M));
	if ( !sp->SfrontM ) {
	  S_free(sp);
	  return 1;
	}
      }
    }
  }
  /*
   *    extend size of V/Vf/VfrontN+M
   */
  if ( sp->flags&S_UVTABLE ) {
    if ( (sp->flags&S_FLOAT)==0 ) {
      sp->V = realloc(sp->V,sizeof(sp->V[0])*N);
      if ( !sp->V ) {
	S_free(sp);
	return 1;
      }
    } else {
      sp->Vf = realloc(sp->Vf,sizeof(sp->Vf[0])*N);
      if ( !sp->Vf ) {
	S_free(sp);
	return 1;
      }
      if ( N-M > sp->usedN-sp->usedM ) {
	/*
	 *   stores *last* values, so cannot shrink or loose a few
	 */
	sp->VfrontM = realloc(sp->VfrontM,sizeof(sp->VfrontM[0])*(N-M+1));
	if ( !sp->VfrontM ) {
	  S_free(sp);
	  return 1;
	}
      }
    }
  }
  /*
   *    now create new vectors
   */
  for (n=sp->usedN+1; n<=N; n++) {
    if ( sp->flags&S_STABLE ) {
      if ( (sp->flags&S_FLOAT)==0 ) {
	sp->S[n-3] = malloc(sizeof(sp->S[0][0])*(M-1));
	if ( !sp->S[n-3] ) {
	  S_free(sp);
	  return 1;
	}
      } else {
	sp->Sf[n-3] = malloc(sizeof(sp->Sf[0][0])*(M-1));
	if ( !sp->Sf[n-3] ) {
	  S_free(sp);
	  return 1;
	}
      }
    }
    if ( sp->flags&S_UVTABLE ) {
      if ( (sp->flags&S_FLOAT)==0 ) {
	sp->V[n-2] = malloc(sizeof(sp->V[0][0])*(M-1));
	if ( !sp->V[n-2] ) {
	  S_free(sp);
	  return 1;
	}
      } else {
	sp->Vf[n-2] = malloc(sizeof(sp->Vf[0][0])*(M-1));
	if ( !sp->Vf[n-2] ) {
	  S_free(sp);
	  return 1;
	}
      }
    }
  }
  {
    int oldN, oldM;
    oldN = sp->usedN;
    oldM = sp->usedM;
    sp->usedN = N;
    sp->usedM = M;
    return S_remake_part(sp,sp->a,oldN,oldM);
  }
}

/*
 *   using cache if possible,
 *   has its own maximum memory allowed to be
 *   greater than usedN
 */
double S_S1(stable_t *sp, unsigned n) {
  if ( n==0 )
    return -INFINITY;
  if ( !sp->S1 )
    return -INFINITY;
  if ( n>sp->usedN ) {
    if ( n>sp->usedN1 ) {
      /*
       *   possibly extend memory and initialise
       */
      int i;
      if ( n>sp->maxN ) {
	return -INFINITY;
      }
      if ( n<sp->usedN1*1.1 )
	n = sp->usedN1*1.1;
      if ( n<sp->usedN1+50 )
	n = sp->usedN1 + 50;
      if ( n>sp->maxN )
	n = sp->maxN;
      sp->S1 = realloc(sp->S1, sizeof(sp->S1[0])*n);
      if ( !sp->S1 ) {
	return -INFINITY;
      }
      for (i=sp->usedN1; i<n; i++) 
	sp->S1[i] = 0;
      sp->usedN1 = n; 
    }
    if ( sp->S1[n-1]==0 ) {
      if ( sp->S1[n-2]==0 )
	sp->S1[n-1] = lgamma(n-sp->a) - sp->lga;
      else
	sp->S1[n-1] = sp->S1[n-2] + log(n-1-sp->a);
    }
  }
  return sp->S1[n-1];
}

double S_U(stable_t *sp, unsigned n, unsigned m) {
  if ( m==1 )
    return n - sp->a;
  assert(m>1);
  return n - m*sp->a + 1/S_V(sp,n,m);
}

double S_UV(stable_t *sp, unsigned n, unsigned m) {
  double SV;
  if ( m==1 )
    return -INFINITY;
  if ( m==n+1 )
	return 1;
  SV = S_V(sp,n,m);
  return (n - m*sp->a)*SV + 1.0;
}

 
double S_V(stable_t *sp, unsigned n, unsigned m) {
  if ( (sp->flags & S_UVTABLE)==0 )
    return 0;
   if ( m>sp->usedM || n>sp->usedN ) {
    if ( n>sp->maxN || m>sp->maxM  ) {
      if ( (sp->flags & S_QUITONBOUND) ) {
	assert(n>sp->maxN || m>sp->maxM);
	if ( sp->tag )
	  yaps_quit("S_V(%u,%u,%lf) tagged '%s' hit bounds\n",n,m,sp->a,
		    sp->tag);
	else 
	  yaps_quit("S_V(%u,%u,%lf) hit bounds\n",n,m,sp->a);
      } else
	return 0;
    }
    if ( S_extend(sp,n,m) ) 
      yaps_quit("S_extend() out of memory\n");
  }
  assert(m>=2);
  if ( n<m ) return 0;
  if ( sp->flags & S_FLOAT ) {
    assert(sp->Vf);
    assert(sp->Vf[n-2]);
    return sp->Vf[n-2][m-2];
  }
  assert(sp->V);
  assert(sp->V[n-2]);
  return sp->V[n-2][m-2];
}

double S_S(stable_t *sp, unsigned N, unsigned T) {
  if ( (sp->flags & S_STABLE)==0 )
    return -INFINITY;
  if ( T>sp->usedM || N>sp->usedN ) {
    if ( N>sp->maxN || T>sp->maxM  ) {
      if ( (sp->flags & S_QUITONBOUND) )
	if ( sp->tag )
	  yaps_quit("S_S(%u,%u,%lf) tagged '%s' hit bounds\n",N,T,sp->a,
		    sp->tag);
	else 
	  yaps_quit("S_S(%u,%u,%lf) hit bounds\n",N,T,sp->a);
      else
	return -INFINITY;
    }
    if ( S_extend(sp,N,T) )
      yaps_quit("S_extend() out of memory\n");
  }
  if ( N==T )
    return 0;
  if ( T==1 )
    return S_S1(sp, N);
  if ( N<T || T==0 )
    return -INFINITY;
  if ( sp->flags&S_FLOAT ) {
    assert(sp->Sf);
    assert(sp->Sf[N-3]);
    return sp->Sf[N-3][T-2];
  } 
  assert(sp->S);
  assert(sp->S[N-3]);
  return sp->S[N-3][T-2];
}

/*
 *    only frees allocated memory, can be called
 *    during a failed make too
 */
void S_free(stable_t *sp) {
  if ( !sp )
    return;
  if ( sp->S1 )  free(sp->S1);
  if ( sp->SfrontN )  free(sp->SfrontN);
  if ( sp->SfrontM )  free(sp->SfrontM);
  if ( sp->VfrontN )  free(sp->VfrontN);
  if ( sp->VfrontM )  free(sp->VfrontM);
  if ( sp->S ) {
    int n;
    for (n=sp->startM-2; n<=sp->usedN-3; n++)
      if ( sp->S[n] ) free(sp->S[n]);
    free(sp->S[0]);
    free(sp->S);
  }
  if ( sp->Sf ) {
    int n;
    for (n=sp->startM-2; n<=sp->usedN-3; n++)
      if ( sp->Sf[n] ) free(sp->Sf[n]);
    free(sp->Sf[0]);
    free(sp->Sf);
  }
  if ( sp->V ) {
    int n;
    for (n=sp->startM-1; n<=sp->usedN-2; n++)
      if ( sp->V[n] ) free(sp->V[n]);
    free(sp->V[0]);
    free(sp->V);
  }
  if ( sp->Vf ) {
    int n;
    for (n=sp->startM-1; n<=sp->usedN-2; n++)
      if ( sp->Vf[n] ) free(sp->Vf[n]);
    free(sp->Vf[0]);
    free(sp->Vf);
  }
  free(sp);
}

void S_report(stable_t *sp, FILE *fp) {
  if ( sp->tag ) 
    fprintf(fp, "Stable '%s': ", sp->tag);
  else
    fprintf(fp, "Stable: ");
  fprintf(fp, "a=%lf, N=%u/%u, M=%u/%u, %s%s %s",
	  sp->a, sp->usedN, sp->maxN, sp->usedM, sp->maxM, 
	  (sp->flags&S_STABLE)?"+S":"", 
	  (sp->flags&S_UVTABLE)?"+U/V":"", 
	  (sp->flags&S_FLOAT)?"float":"double");
#ifdef MEMALLOCED
  fprintf(fp, " mem=%uk\n", sp->memalloced/1024);
#endif
  fprintf(fp, "\n");
}
