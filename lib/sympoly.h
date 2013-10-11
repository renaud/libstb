/*
 * Elementary symmetric polynomials
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
 *  See http://en.wikipedia.org/wiki/Newton%27s_identities
 *  for math details of elementary symmetric polynomials,
 *  but we don't use their recursions because they are unstable.
 *  
 *  Based on an alternative recursion (on K):
 *      e_H(x_0,...,x_{K-1}) = 
 *            1_{H<K} && e_H(x_0,...,x_{K-2}) + x_{K-1}e_{H-1}(x_0,...,x_{K-2})
 *  to build values in situ in a vector or matrix.
 *  this recursion is more stable and allows a simple trick to
 *  deal with overflow.  If many x_k values >1, can get overflow.
 *  So instead we collect these together and run the recursion on:
 *  f_H(x_0,...,x_{K-1}) == \prod_{k=0}{^{K-1} (1/x_k)^{x_k>1}
 *                               * e_H(x_0,...,x_{K-1})
 */

#ifndef __SYMPOLY_H
#define __SYMPOLY_H

#include <stdint.h>
#include "rng.h"

/*
 *   for optional statically declared arrays, but these are just
 *   buffers that are worked around using malloc(), so no need
 *   to change ....
 */
#define SYMPOLY_MAX 10

/*
`*  Compute elementary symmetric polynomials
 *     res[0] = 1
 *     res[1] = e_1(x_0,...,x_{K-1})
 *     res[2] = e_2(x_0,...,x_{K-1})
 *     ...
 *     res[K] = e_K(...)
 *
 *  is O(K^2) 
 *
 *  if values get too big, log overflow will appear in (*overflow),
 *    reconstruct values as
 *       res[0], exp(*overflow)*(res[1],...,res[K])
 *  otherwise  overflow=0
 *
 *  return non-zero on error (e.g.,  K too big)
 */
int sympoly(int K, double *val, double *res, double *overflow);

/*
 *    sample exactly H of the K features proportional to
 *    occurrences in the elem.sym.poly.
 *
 *    is O(HK) for H<=K
 *
 *    results returned as a bit-vector, 
 *           i.e.,  max. for H is 31!
 *
 *    return 0 if OK, non-zero on error
 */
uint32_t sympoly_sample(int K, int H, double *val, rngp_t rng);
 
#endif
