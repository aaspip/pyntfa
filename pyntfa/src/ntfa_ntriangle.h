/* Non-stationary triangle smoothing */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
// #include "_bool.h"
// #include "alloc.h"

#include <stdbool.h>
#include "ntfa_alloc.h"

#ifndef _ntriangle_h

typedef struct tf_Ntriangle *tf_ntriangle;
/* abstract data type */
/*^*/

#endif

struct tf_Ntriangle {
    float *tmp;
    int np, nb, nx;
};

// static void fold (int o, int d, int nx, int nb, int np, 
// 		  const float *x, float* tmp);
// static void fold2 (int o, int d, int nx, int nb, int np, 
// 		   float *x, const float* tmp);
// static void doubint (int nx, float *x, bool der);
// static void triple (int o, int d, int nx, int nb, 
// 		    const float* t, const float* s, float* x, const float* tmp);
// static void triple2 (int o, int d, int nx, int nb, 
// 		     const float* t, const float* s, const float* x, float* tmp);
// static void double1 (int o, int d, int nx, int nb, 
//         const float* t, const float* s, float* x, const float* tmp);

tf_ntriangle tf_ntriangle_init (int nbox /* maximum triangle length */, 
			     int ndat /* data length */);
/*< initialize >*/

void tf_nsmooth (tf_ntriangle tr   /* smoothing object */, 
		 int o, int d   /* sampling */, 
		 bool der       /* derivative flag */, 
		 const float *t /* triangle lengths */, 
		 const float *s /* triangle shifts */,
		 float *x       /* data (smoothed in place) */);
/*< smooth >*/

void tf_nsmooth2 (tf_ntriangle tr   /* smoothing object */, 
		  int o, int d   /* sampling */, 
		  bool der       /* derivative flag */, 
		  const float *t /* triangle lengths */,
		  const float *s /* triangle shifts */,
		  float *x       /* data (smoothed in place) */);
/*< alternative smooth >*/

void tf_ndsmooth (tf_ntriangle tr /* smoothing derivative object */, 
		  int o, int d /* sampling. o: starting index, d: stride in samples for 1/2/3rd dimension; all refer to a correct index in a 1D vector  */, 
		  bool der     /* derivative flag */, 
		  const float *t /* triangle lengths */, 
		  const float *s /* triangle shifts */,
		  float *x     /* data (smoothed in place) */);
/*< smooth derivative >*/

void  tf_ntriangle_close(tf_ntriangle tr);
/*< free allocated storage >*/

