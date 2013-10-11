/*
 * caching difference of lgamma and digamma
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
 *     
 */

#ifndef __LGAMMA_H
#define __LGAMMA_H

/*
 *   caches are self allocating
 */
#define GCACHE 100
struct gcache_s {
  double par;
  double lgpar;
  double cache[GCACHE];
} ;

void gcache_init(struct gcache_s *lpg, double p);
double gcache_value(struct gcache_s *lpg, int j);
void pcache_init(struct gcache_s *lpg, double p);
double pcache_value(struct gcache_s *lpg, int j);
void qcache_init(struct gcache_s *lpg, double p);
double qcache_value(struct gcache_s *lpg, int j);

double gammadiff(int N, double alpha, double lga);
double psidiff(int N, double alpha, double pa);

#endif

