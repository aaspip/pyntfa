

#include <stdbool.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include "ntfa_alloc.h"
// #include "ntfa_conjgrad.h"
// #include "ntfa_ntriangle.h"
// #include "ntfa_ntrianglen.h"
// #include "ntfa_decart.h"

void tf_ntrianglen_init (int ndim    /* number of dimensions */, 
		      int *nbox   /* triangle radius [ndim] */, 
		      int *ndat   /* data dimensions [ndim] */,
		      float **len /* triangle lengths [ndim][nd] */,
                      float **sft /* triangle shifts [ndim][nd] */,
		      int repeat  /* repeated smoothing */);
/*< initialize >*/

void tf_ntrianglen_lop (bool adj, bool add, int nx, int ny, float* x, float* y);
/*< linear operator >*/

void tf_ntrianglen_close(void);
/*< free allocated storage >*/
