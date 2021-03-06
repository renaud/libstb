
ABOUT
=====

This source directory contains library routines that provide alternative ways of computing generalised second order Stirling numbers used in working with Pitman-Yor and Dirichlet processes (PYP and DP). Included is library routines for posterior sampling on the discount and concentration parameters of the PYP/DP, and some simple demo scripts.

System used for Pitman-Yor processing with topic models, to allow easy scaling to gigabytes of text. Tested on a few versions on Ubuntu Linux and MacOSX. Requires the GSL (GNU Scientific Library) for polygamma functions.

Matching paper with some of the theory at (http://arxiv.org/abs/1007.0296), and a series of papers by Buntine and students. 


CONTENTS
========
This directory contains a library `lib/` and some testing programs `test/`.

The library provides alternative ways of computing generalised
second order Stirling numbers used in working with Pitman-Yor
and Dirichlet processes (PYP and DP).  The relevant theory appears in:
    http://arxiv.org/abs/1007.0296v2
and some additional theory is given in "doc/alpha.pdf".

The library functions are described in the header files:
    stable.h   - tables for log of Stirling numbers
    		  and for ratio of Stirling numbers
    psample.h  - routines for sampling discount (a) and concentration (b)
    sapprox.h  - approximate calcs from difference approximation
    arms.h     - the ARM code from Gilks, used by samplers
    sympoly.h  - evaluation and sampling of elementary symmetric polynomials
		 (useful for some more complex models)

In addition the following files give links to key external 
libraries that you can change if you want to alter dependencies:
    rng.h digamma.h  
For testing purposes, the following is kept:   yaps.h
An earlier full featured version is also here,
roughly documented.  Used for really big tables:  sxtable.h
Note am slowly integrating the functionality of this into the
better engineered stable.h.

Including is testers, "list.c" and "demo.c", 
an example use as well as a worked
example of posterior sampling on the discount and concentration
parameters of the PYP/DP, in psample.h.

The tester has all sorts of dependencies, but the library
functions are self-contained apart from requiring 
digamma() and polygamma() functions 
(currently got from GNU Scientific Library, GSL)
and random number generators.
These dependencies are isolated in "lib/digamma.h" and
"lib/rng.h" so redefine these to change these around.

Building
=======
Tested on a few versions on Ubuntu Linux and an old MacOSX.
If "lib/digamma.h" and "lib/rng.h" are unchanged you will need
to have GSL installed.
Testing also requires the GSL for randon number functions.
    NB.  since its good to have control of seeds 
Note to compile faster versions, modify your compile time flags
accordingly,  e.g.,  -O5 -DNDEBUG
MacOSX requires fiddling with paths for GSL libraries.
   cd lib
   #  calling plain make creates the library, and leaves it at ./libstb.a
   make

Optional
========
   cd test
   #  - a simple value lister, lets you see the values
   make list
   #  - see what it does
   ./list
   #  - create the "demo" executable with all sorts of bells and whistles
   #    for evaluating the Stirling numbers
   #  - to see the options, execute:
   #        ./demo -h
   make demo
   #  - check the simple options
   ./demo -h
   #  - run with fitting the concentration
   ./demo -H5 -b 10

Installing
=========
Copy the library "lib/libstb.a" and the header files 
"stable.h", "sapprox.h", to wherever is useful for you.
If using the parameter samplers, samplea() and sampleb(), then
you also need "psample.h", "arms.h", "rng.h" and "lgamma.h".

Detail
======
To see how to use the library in more detail, see the examples
in "list.c" and "demo.c".  This also shows how to do posterior sampling
on the discount and concentration parameters of the PYP/DP.

In addition, "precision_test.c" is a stand alone program to
compare the two Stirling number recursions when done in
float versus double.  This is intended to evaluate how the
finite precision affects the calculations, and demonstrates
the "ratio of Stirling number" recursion is far more
accurate.

The adaptive rejection sampling code from:
   http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html
is included as "arms.h" and "arms.c".
This is optional since the slice sampler is the default.
Note naive slice sampling *doesn't* work when sampling the concentration
parameter ... you need to warm up the sampler with a bit of fixed
point optimisation.

See the relevant "#define" in "lib/psample.h".

LICENSE
=======
The code in this directory is subject to the Mozilla Public License Version 1.1 
(the "License"); you may not use files in this directory except in
compliance with the License. You may obtain a copy of the License at
http://www.mozilla.org/MPL/

Software distributed under the License is distributed on an "AS IS"
basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
License for the specific language governing rights and limitations
under the License.

Authors:
     Wray Buntine (wray.buntine@nicta.com.au)
     Lan Du (lan.du@nicta.com.au)
Contact: 
     Wray Buntine (wray.buntine@nicta.com.au)


