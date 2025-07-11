
#ifndef _tf_triangle_h
#define _tf_triangle_h


// #include "_bool.h"


typedef struct tf_Triangle *tf_triangle;
/* abstract data type */


tf_triangle tf_triangle_init (int  nbox /* triangle length */, 
			      int  ndat /* data length */,
                              bool box  /* if box instead of triangle */);
/*< initialize >*/


void tf_smooth (tf_triangle tr  /* smoothing object */, 
		int o, int d    /* trace sampling */, 
		bool der        /* if derivative */, 
		float *x        /* data (smoothed in place) */);
/*< apply triangle smoothing >*/


void tf_dsmooth (tf_triangle tr  /* smoothing object */, 
		int o, int d    /* trace sampling */, 
		bool der        /* if derivative */, 
		float *x        /* data (smoothed in place) */);
/*< apply triangle smoothing >*/


void tf_smooth2 (tf_triangle tr  /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */,
		 float *x        /* data (smoothed in place) */);
/*< apply adjoint triangle smoothing >*/


void tf_dsmooth2 (tf_triangle tr  /* smoothing object */, 
		 int o, int d    /* trace sampling */, 
		 bool der        /* if derivative */,
		 float *x        /* data (smoothed in place) */);
/*< apply adjoint triangle smoothing >*/


void  tf_triangle_close(tf_triangle tr);
/*< free allocated storage >*/

#endif
