#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "ntfa_divnnsc.h"

#define SF_MAX(a,b) ((a) < (b) ? (b) : (a))
#define SF_MIN(a,b) ((a) < (b) ? (a) : (b))
#define SF_MAX_DIM 9
#define SF_PI (3.14159265358979323846264338328)

#include "ntfa_alloc.h"

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

/*from user/chenyk/vecoper.c */
static PyObject *ntft(PyObject *self, PyObject *args){
    
	/**initialize data input**/
    int nd;
    
    PyObject *f1=NULL;
    PyObject *arrf1=NULL;

	int ndata;	/*integer parameter*/
	float fpar; /*float parameter*/
    int ndim, i;
    float *data;
    
	PyArg_ParseTuple(args, "Oif", &f1, &ndata, &fpar);

    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    
    npy_intp *sp=PyArray_SHAPE(arrf1);

	data  = (float*)malloc(ndata * sizeof(float));
	
    if (*sp != ndata)
    {
    	printf("Dimension mismatch, N_input = %d, N_data = %d\n", *sp, ndata);
    	return NULL;
    }
    
    /*reading data*/
    for (i=0; i<ndata; i++)
    {
        data[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }

	/*sub-function goes here*/
	
	
    /*Below is the output part*/
    PyArrayObject *vecout;
	npy_intp dims[2];
	dims[0]=ndata;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = data[i];
	
	return PyArray_Return(vecout);
	
}


static PyObject *tf1d(PyObject *self, PyObject *args){
    
	/**initialize data input**/
    int nd;
    
    PyObject *f1=NULL;
    PyObject *arrf1=NULL;

	int ndata;	/*integer parameter*/
	float fpar; /*float parameter*/
    int ndim, i;
    float *data;
    
    int niter,verb,rect0,n1;
    float dt,alpha;
    int ifb,inv;
	PyArg_ParseTuple(args, "Oiiiiiiff", &f1,&niter,&n1,&verb,&rect0,&ifb,&inv,&dt,&alpha);
	
	ndata=n1;


    int i1, iw, nt, nw, i2, n2=1, n12, n1w;
    int m[SF_MAX_DIM], *rect;
    float t, d1, w, w0, dw, mean=0.0f;
    float *trace, *kbsc, *mkbsc, *sscc, *mm, *ww;
    
    d1=dt;
	    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
	    nw = nt/2+1;
	    dw = 1./(nt*d1);
	    w0 = 0.;


	for(i2=0; i2 < SF_MAX_DIM; i2 ++) {
	    m[i2] = 1;
	}
	m[0] = n1;

    n1w = n1*nw;
    n12 = 2*n1w;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;


	printf("niter=%d,n1=%d,nw=%d,nt=%d,n12=%d,rect0=%d\n",niter,n1,nw,nt,n12,rect0);
	printf("dw=%g,w0=%g,alpha=%g,dt=%g\n",dw,w0,alpha,dt);

	

    trace = tf_floatalloc(n1);
    kbsc    = tf_floatalloc(n12); /*kbsc is the basis functions*/

    rect = tf_intalloc(2*nw);
    for (iw=0; iw < nw; iw++) {
	rect[iw+nw] = rect[iw] = SF_MAX(1, (int) rect0/(1.0+alpha*iw/nw));
    }

	sscc = tf_floatalloc(n12);
	divnn_sc_init(2*nw, 1, n1, m, rect, kbsc, 
			(bool) (verb && (n2 < 500)));  /* Compared with the non-stationary version, only difference*/
			
	mm = NULL;
	ww = NULL;

    /* sin and cos basis */
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    if (0.==w) { /* zero frequency */
		kbsc[iw*n1+i1] = 0.;
	    } else {
		t = i1*d1;
		kbsc[iw*n1+i1] = sinf(w*t);
	    }
	}
    }
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    if (0.==w) { /* zero frequency */
		kbsc[(iw+nw)*n1+i1] = 0.5;
	    } else {
		t = i1*d1;
		kbsc[(iw+nw)*n1+i1] = cosf(w*t); /*YC 10/20/2024: kbsc is the basis functions*/
	    }
	}
    }
    
    if (NULL != mm || NULL != ww) {
	mkbsc = tf_floatalloc(n12);
	for (i1=0; i1 < n12; i1++) {
	    mkbsc[i1] = kbsc[i1];
	}
    } else {
	mkbsc = NULL;

	mean = 0.;
	for (i1=0; i1 < n12; i1++) {
	    mean += kbsc[i1]*kbsc[i1];
	}
	mean = sqrtf (mean/(n12));
	for (i1=0; i1 < n12; i1++) {
	    kbsc[i1] /= mean;
	}
    }
    
	if (NULL != mm || NULL != ww) {
	    for (i1=0; i1 < n12; i1++) {
		kbsc[i1] = mkbsc[i1];
	    }

	    if (NULL != mm) {
		for (iw=0; iw < 2*nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			kbsc[iw*n1+i1] *= mm[i1];
		    }
		}
	    }

	    if (NULL != ww) {
		for (iw=0; iw < nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			kbsc[iw*n1+i1]      *= ww[iw*n1+i1];
			kbsc[(iw+nw)*n1+i1] *= ww[iw*n1+i1];
		    }
		}
	    }

	    mean = 0.;
	    for (i1=0; i1 < n12; i1++) {
		mean += kbsc[i1]*kbsc[i1];
	    }
	    mean = sqrtf (mean/(n12));
	    for (i1=0; i1 < n12; i1++) {
		kbsc[i1] /= mean;
	    }
	}
	    
	if(!inv)
	{
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    npy_intp *sp=PyArray_SHAPE(arrf1);

	data  = (float*)malloc(ndata * sizeof(float));
	
    if (*sp != ndata)
    {
    	printf("Dimension mismatch, N_input = %d, N_data = %d\n", *sp, ndata);
    	return NULL;
    }
    
    /*reading data*/
    for (i=0; i<ndata; i++)
    {
        trace[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }
	printf("ndata=%d\n",ndata);
	    
	    if (NULL != mm) {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] *= mm[i1];
		}
	    }
	    
	    for(i1=0; i1 < n1; i1++) {
		trace[i1] /= mean;
	    }
	    divnn_sc (trace,sscc,niter); /*YC 10/20/2024: sscc is the NTF coefficients*/
	}else{
	/*This part is to reconstruct the data given the basis functions and their weights (i.e., TF spectrum)*/
	
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    npy_intp *sp=PyArray_SHAPE(arrf1);
    	
    for (i=0; i<n1w*2; i++)
    {
        sscc[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }

	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = 0.;
	    }
	    
	    for (iw=0; iw < nw; iw++) {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] += sscc[(iw+nw)*n1+i1]*kbsc[(iw+nw)*n1+i1]
			*mean+sscc[iw*n1+i1]*kbsc[iw*n1+i1]*mean;
		    if (NULL != mm) trace[i1] *= mm[i1];
		}
	    }
	}
	/*sub-function goes here*/
	
    /*Below is the output part*/
    PyArrayObject *vecout;
    npy_intp dims[2];
    
    if(!inv)
    {

	if(ifb) /*if output the basis functions, e.g., Fourier bases in this case*/
	{
	dims[0]=ndata*nw*4+3;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = sscc[i];
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i+ndata*nw*2))) = kbsc[i];
		
	(*((float*)PyArray_GETPTR1(vecout,0+ndata*nw*4))) = w0;
	(*((float*)PyArray_GETPTR1(vecout,1+ndata*nw*4))) = dw;
	(*((float*)PyArray_GETPTR1(vecout,2+ndata*nw*4))) = nw;
	}else{
	dims[0]=ndata*nw*2+3;dims[1]=1;
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = sscc[i];
	printf("w0=%g,dw=%g,nw=%d\n",w0,dw,nw);
	(*((float*)PyArray_GETPTR1(vecout,0+ndata*nw*2))) = w0;
	(*((float*)PyArray_GETPTR1(vecout,1+ndata*nw*2))) = dw;
	(*((float*)PyArray_GETPTR1(vecout,2+ndata*nw*2))) = nw;
	}
	
	free(kbsc);free(sscc);free(trace);free(rect);
	divnn_sc_close();
	
	}else{
	
	dims[0]=n1;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = trace[i];
	
	free(trace);
	}
	
	return PyArray_Return(vecout);
	
}

static PyObject *ntf1d(PyObject *self, PyObject *args){
    
	/**initialize data input**/
    int nd, nd2;
    
    PyObject *f1=NULL;
    PyObject *arrf1=NULL;
    PyObject *f2=NULL;
    PyObject *arrf2=NULL;
    PyObject *f3=NULL;
    PyObject *arrf3=NULL;
    
	int ndata;	/*integer parameter*/
	float fpar; /*float parameter*/
    int ndim, i;
    float *data;
    
    int niter,verb,rect0,n1;
    float dt,alpha;
    int ifb,inv;
	PyArg_ParseTuple(args, "OOOiiiiiiff", &f1,&f2,&f3,&niter,&n1,&verb,&rect0,&ifb,&inv,&dt,&alpha);
	
	ndata=n1;


    int i1, iw, nt, nw, i2, n2=1, n12, n1w;
    int m[SF_MAX_DIM], *rect;
    float t, d1, w, w0, dw, mean=0.0f;
    float *trace, *kbsc, *mkbsc, *sscc, *mm, *ww;
    int box[SF_MAX_DIM];
    
    	d1=dt;
	    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
	    nw = nt/2+1;
	    dw = 1./(nt*d1);
	    w0 = 0.;


	for(i2=0; i2 < SF_MAX_DIM; i2 ++) {
	    m[i2] = 1;
	}
	m[0] = n1;

    n1w = n1*nw;
    n12 = 2*n1w;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

	printf("running ntf1d\n");
	printf("niter=%d,n1=%d,nw=%d,nt=%d,n12=%d,rect0=%d\n",niter,n1,nw,nt,n12,rect0);
	printf("dw=%g,w0=%g,alpha=%g,dt=%g\n",dw,w0,alpha,dt);

    float *rct[2]; /* storing non-stationary smoothing radii */
    
	rct[0] = tf_floatalloc (n1*nw*2);
	rct[1] = tf_floatalloc (n1*nw*2);
	
    float  *sft[2];	/* storing non-stationary shifting size */
	sft[0] = tf_floatalloc (n1*nw*2);
	sft[1] = tf_floatalloc (n1*nw*2);

    trace = tf_floatalloc(n1);
    kbsc    = tf_floatalloc(n12); /*kbsc is the basis functions*/

	for(i=0;i<n1;i++) trace[i]=0;
	for(i=0;i<n12;i++) kbsc[i]=0;

    rect = tf_intalloc(2*nw);
    for (iw=0; iw < nw; iw++) {
	rect[iw+nw] = rect[iw] = SF_MAX(1, (int) rect0/(1.0+alpha*iw/nw));
    }

	sscc = tf_floatalloc(n12);
	
    arrf2 = PyArray_FROM_OTF(f2, NPY_FLOAT, NPY_IN_ARRAY);
    arrf3 = PyArray_FROM_OTF(f3, NPY_FLOAT, NPY_IN_ARRAY);

	for(i=0;i<n1*nw;i++)
	{
	rct[0][i]=*((float*)PyArray_GETPTR1(arrf2,i));
	rct[0][n1*nw+i]=rct[0][i];
	rct[1][i]=*((float*)PyArray_GETPTR1(arrf3,i));
	rct[1][n1*nw+i]=rct[1][i];
	sft[0][i]=0;
	sft[0][n1*nw+i]=0;
	sft[1][i]=0;
	sft[1][n1*nw+i]=0;
	}
	
/*	divnn_sc_init(2*nw, 1, n1, m, rect, kbsc, 
			(bool) (verb && (n2 < 500))); */ /* Compared with the stationary version, only difference*/
	m[0]=n1;
	m[1]=nw;
	m[2]=2;
	
	box[0]=rct[0][0];box[1]=rct[1][0];
	for(i=0;i<n1*nw;i++)
	{
	if(box[0]<rct[0][i]) box[0]=rct[0][i];
	if(box[1]<rct[1][i]) box[1]=rct[1][i];
	}
	box[2]=1;
	printf("box[0]=%d,box[1]=%d,box[2]=%d\n",box[0],box[1],box[2]);
	
	divnn_sc_init2(2*nw, 3, n1, m, box, rct,sft, kbsc, 
			(bool) (verb && (n2 < 500))); 
	
	mm = NULL;
	ww = NULL;

    /* sin and cos basis */
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    if (0.==w) { /* zero frequency */
		kbsc[iw*n1+i1] = 0.;
	    } else {
		t = i1*d1;
		kbsc[iw*n1+i1] = sinf(w*t);
	    }
	}
    }
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    if (0.==w) { /* zero frequency */
		kbsc[(iw+nw)*n1+i1] = 0.5;
	    } else {
		t = i1*d1;
		kbsc[(iw+nw)*n1+i1] = cosf(w*t); /*YC 10/20/2024: kbsc is the basis functions*/
	    }
	}
    }
    
    if (NULL != mm || NULL != ww) {
	mkbsc = tf_floatalloc(n12);
	for (i1=0; i1 < n12; i1++) {
	    mkbsc[i1] = kbsc[i1];
	}
    } else {
	mkbsc = NULL;

	mean = 0.;
	for (i1=0; i1 < n12; i1++) {
	    mean += kbsc[i1]*kbsc[i1];
	}
	mean = sqrtf (mean/(n12));
	for (i1=0; i1 < n12; i1++) {
	    kbsc[i1] /= mean;
	}
    }
    
	if (NULL != mm || NULL != ww) {
	    for (i1=0; i1 < n12; i1++) {
		kbsc[i1] = mkbsc[i1];
	    }

	    if (NULL != mm) {
		for (iw=0; iw < 2*nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			kbsc[iw*n1+i1] *= mm[i1];
		    }
		}
	    }

	    if (NULL != ww) {
		for (iw=0; iw < nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			kbsc[iw*n1+i1]      *= ww[iw*n1+i1];
			kbsc[(iw+nw)*n1+i1] *= ww[iw*n1+i1];
		    }
		}
	    }

	    mean = 0.;
	    for (i1=0; i1 < n12; i1++) {
		mean += kbsc[i1]*kbsc[i1];
	    }
	    mean = sqrtf (mean/(n12));
	    for (i1=0; i1 < n12; i1++) {
		kbsc[i1] /= mean;
	    }
	}
	    
	if(!inv)
	{
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    npy_intp *sp=PyArray_SHAPE(arrf1);

	data  = (float*)malloc(ndata * sizeof(float));
	
    if (*sp != ndata)
    {
    	printf("Dimension mismatch, N_input = %d, N_data = %d\n", *sp, ndata);
    	return NULL;
    }
    
    /*reading data*/
    for (i=0; i<ndata; i++)
    {
        trace[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }
	
	    if (NULL != mm) {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] *= mm[i1];
		}
	    }
	    
	    for(i1=0; i1 < n1; i1++) {
		trace[i1] /= mean;
	    }
	    /*divnn_sc (trace,sscc,niter);*/ /*YC 10/20/2024: sscc is the NTF coefficients*/
	    divnn_sc2 (trace,sscc,niter); /*Compared with divnn_sc, no big difference*/

	}else{
	/*This part is to reconstruct the data given the basis functions and their weights (i.e., TF spectrum)*/
	
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    npy_intp *sp=PyArray_SHAPE(arrf1);
    	
    for (i=0; i<n1w*2; i++)
    {
        sscc[i]=*((float*)PyArray_GETPTR1(arrf1,i));
    }

	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = 0.;
	    }
	    
	    for (iw=0; iw < nw; iw++) {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] += sscc[(iw+nw)*n1+i1]*kbsc[(iw+nw)*n1+i1]
			*mean+sscc[iw*n1+i1]*kbsc[iw*n1+i1]*mean;
		    if (NULL != mm) trace[i1] *= mm[i1];
		}
	    }
	}
	/*sub-function goes here*/
	
    /*Below is the output part*/
    PyArrayObject *vecout;
    npy_intp dims[2];
    
    if(!inv)
    {

	if(ifb) /*if output the basis functions, e.g., Fourier bases in this case*/
	{
	dims[0]=ndata*nw*4+3;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = sscc[i];
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i+ndata*nw*2))) = kbsc[i];
		
	(*((float*)PyArray_GETPTR1(vecout,0+ndata*nw*4))) = w0;
	(*((float*)PyArray_GETPTR1(vecout,1+ndata*nw*4))) = dw;
	(*((float*)PyArray_GETPTR1(vecout,2+ndata*nw*4))) = nw;
	}else{
	dims[0]=ndata*nw*2+3;dims[1]=1;
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<ndata*nw*2;i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = sscc[i];
	printf("w0=%g,dw=%g,nw=%d\n",w0,dw,nw);
	(*((float*)PyArray_GETPTR1(vecout,0+ndata*nw*2))) = w0;
	(*((float*)PyArray_GETPTR1(vecout,1+ndata*nw*2))) = dw;
	(*((float*)PyArray_GETPTR1(vecout,2+ndata*nw*2))) = nw;
	}

	free(rct[0]);free(rct[1]);
	free(sft[0]);free(sft[1]);
	free(kbsc);free(sscc);free(trace);free(rect);
	divnn_sc_close2();
	
	}else{
	
	dims[0]=n1;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = trace[i];
		
	free(trace);
	}

	
	
	return PyArray_Return(vecout);

}

/*documentation for each functions.*/
static char ntfacfun_document[] = "Document stuff for this C module...";

/*defining our functions like below:
  function_name, function, METH_VARARGS flag, function documents*/
static PyMethodDef functions[] = {
  {"ntft", ntft, METH_VARARGS, ntfacfun_document},
  {"tf1d", tf1d, METH_VARARGS, ntfacfun_document},
  {"ntf1d", ntf1d, METH_VARARGS, ntfacfun_document},
  {NULL, NULL, 0, NULL}
};

/*initializing our module informations and settings in this structure
for more informations, check head part of this file. there are some important links out there.*/
static struct PyModuleDef ntfacfunModule = {
  PyModuleDef_HEAD_INIT, /*head informations for Python C API. It is needed to be first member in this struct !!*/
  "ntfacfun",  /*module name*/
  NULL, /*means that the module does not support sub-interpreters, because it has global state.*/
  -1,
  functions  /*our functions list*/
};

/*runs while initializing and calls module creation function.*/
PyMODINIT_FUNC PyInit_ntfacfun(void){
  
    PyObject *module = PyModule_Create(&ntfacfunModule);
    import_array();
    return module;
}
