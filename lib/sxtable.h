/*
 * Stirling Number table handling for Pitman-Yor processing
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
 *   Details of the algorithms are in various papers by
 *   Buntine, Du and Hutter.
 *
 *   This is also intended to support Adaptive rejection sampling
 *   of the parameter to the Stirling numbers which means
 *   derivatives and two values (with two full tables) can be dealt with.
 *
 */


/*
 *  access a Stirling Number safely, but gives the log(),
 *  growing the table if needed, and
 *  recursively recomputing it if it is stored sparsely
 */
double S_safe(int N, int M);
/*
 *  above but return derivative of log S w.r.t. a
 */
double S_safe_da(int N, int M);

/*
 *  fills S table with maxN and maxM values;
 *   -- calling S_safe() will check bounds and extend table
 *      automatically if needed
 *   -- set fix!=0 disallows automatic extension of the table
 *   -- setting skip>1 means extensions done will only store
 *      every skip-th of the N rows, hence making the matrix
 *      1/skip-th smaller, but now S evaluation requires
 *      2^skip evaluations since it must be done recursively
 *   -- set divin!=0 if you want to evaluate derivatives too
 *      NB.  require's considerable more memory and complexity
 *           so don't if you don't want it
 */
void S_make(int maxN, int maxM, float a, int skip, int fix, int divin);

/*
 *   return current a value
 */
float S_alpha(void);
/*
 *  does the bounds check and extends if needed, called inside
 *  S_safe() too
 */
void S_check(int N, int M);
void S_report(void);
void S_free(void);
/*
 *  efficient buffering of tables if up to 2 additional
 *  tables should be saved
 */
void S_rebuild(float apar, int primary);
/*
 *  set primary to non-zero if this uses the
 *  preferable value of "a" to keep
 */
void S_save(int I, int primary);
