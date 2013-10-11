/*
 * Approximate S computation
 * Copyright (C) 2009 Wray Buntine 
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
 *  Exact/Approx calculations for S, for theory, see:
 *      http://arxiv.org/abs/1007.0296
 *
 *  Returns -HUGE_VAL if any trouble.
 *         e.g.,  m>4
 *          
 */

/*
 *    approximate calculation for M<=4, exact when a==0
 */
double S_approx(int n, int m, float a);

/*
 *    approx calculation for M<=4, derivative w.r.t. a
 */
double S_approx_da(int n, int m, float a) ;
