/* Non-stationary triangle smoothing */

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

