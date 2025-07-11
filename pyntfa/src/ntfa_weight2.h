
#ifndef _tf_weight2_h
#define _tf_weight2_h


#include <stdbool.h>


void tf_weight2_init(int nw1   /* number of components */, 
		     int n     /* model size */, 
		     float *ww /* weight [nw*n] */);
/*< initialize >*/


void tf_weight2_close(void);
/*< free allocated storage >*/


void tf_weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);
/*< linear operator >*/

#endif

