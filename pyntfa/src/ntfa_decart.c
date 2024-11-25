/* Conversion between line and Cartesian coordinates of a vector. */

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
/*^*/

#include "ntfa_decart.h"

void tf_line2cart(int dim         /* number of dimensions */, 
		  const int* nn /* box size [dim] */, 
		  int i         /* line coordinate */, 
		  int* ii       /* cartesian coordinates [dim] */)
/*< Convert line to Cartesian >*/
{
    int axis;
 
    for (axis = 0; axis < dim; axis++) {
	ii[axis] = i%nn[axis];
	i /= nn[axis];
    }
}

int tf_cart2line(int dim         /* number of dimensions */, 
		 const int* nn /* box size [dim] */, 
		 const int* ii /* cartesian coordinates [dim] */) 
/*< Convert Cartesian to line >*/
{
    int i, axis;

    if (dim < 1) return 0;

    i = ii[dim-1];
    for (axis = dim-2; axis >= 0; axis--) {
	i = i*nn[axis] + ii[axis];
    }
    return i;
}

int tf_first_index (int i          /* dimension [0...dim-1] */, 
		    int j        /* line coordinate */, 
		    int dim        /* number of dimensions */, 
		    const int *n /* box size [dim] */, 
		    const int *s /* step [dim] */)
/*< Find first index for multidimensional transforms >*/
{
    int i0, n123, ii;
    int k;

    n123 = 1;
    i0 = 0;
    for (k=0; k < dim; k++) {
	if (k == i) continue;
	ii = (j/n123)%n[k]; /* to cartesian */
	n123 *= n[k];	
	i0 += ii*s[k];      /* back to line */
    }

    return i0;
}

void tf_large_line2cart(int dim         /* number of dimensions */, 
			const off_t* nn /* box size [dim] */, 
			off_t i         /* line coordinate */, 
			off_t* ii       /* cartesian coordinates [dim] */)
/*< Convert line to Cartesian >*/
{
    int axis;
 
    for (axis = 0; axis < dim; axis++) {
	ii[axis] = i%nn[axis];
	i /= nn[axis];
    }
}

off_t tf_large_cart2line(int dim         /* number of dimensions */, 
			 const off_t* nn /* box size [dim] */, 
			 const off_t* ii /* cartesian coordinates [dim] */) 
/*< Convert Cartesian to line >*/
{
    off_t i;
    int  axis;

    if (dim < 1) return 0;

    i = ii[dim-1];
    for (axis = dim-2; axis >= 0; axis--) {
	i = i*nn[axis] + ii[axis];
    }
    return i;
}

off_t tf_large_first_index (int i          /* dimension [0...dim-1] */, 
			    off_t j        /* line coordinate */, 
			    int dim        /* number of dimensions */, 
			    const off_t *n /* box size [dim] */, 
			    const off_t *s /* step [dim] */)
/*< Find first index for multidimensional transforms >*/
{
    off_t i0, n123, ii;
    int k;

    n123 = 1;
    i0 = 0;
    for (k=0; k < dim; k++) {
	if (k == i) continue;
	ii = (j/n123)%n[k]; /* to cartesian */
	n123 *= n[k];	
	i0 += ii*s[k];      /* back to line */
    }

    return i0;
}


