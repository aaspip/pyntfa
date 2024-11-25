from ntfacfun import *

def tf(din):
	"""
	Time-frequency transform
	
	By Yangkang Chen
	Nov 6, 2022
	"""
	import numpy as np
	
	din=np.float32(din)
	ndata=din.size;
	dout=ntft(din,ndata,3.0)
	
	return dout
	
def ntfa1d(din,dt=0.004,niter=100,rect=10,ifb=0,ifn=0,inv=0,verb=0,alpha=0,rect1=None,rect2=None):
	"""
	Non-stationary time-frequency transform on 1D trace
	
	INPUT
	din: input trace (1D, n1/nt)
	ifb: if output basis functions
	ifn: if non-stationary smoothing for model
	inv: flag of inverse transform (0: forward; 1: inverse; default: 0)
	
	OUTPUT
	dout: output Time-frequency spectrum (first axis: Time; second axis: Frequency; third axis: imag/real)
	NOTE: 1) dout is of size n1*nw*2 (e.g., dout.reshape([n1,nw,2],order='F'))
		  2) The first component in 3rd axis is Imaginary; and the second component is Real (Remember it when constructing a complex number)	
	
	EXAMPLE
	demos/test_pyntfa_syn1d.py
	demos/test_pyntfa_syn1d2.py
	
	HISTORY
	Original version by Yangkang Chen, Nov 6, 2022
	Modified by Yangkang Chen, Oct 20, 2024
	
	"""
	import numpy as np
	
	din=np.float32(din)
	ndata=din.size;
	
# 	if din.ndim==1:	#for 2D problems
# 		n1=din.size;
# 	elif din.ndim==2:
# 		[n1,n2]=din.shape;
	n1=din.shape[0];
	print(n1)
	
	if ifn:
		rect1=np.float32(rect1);rect2=np.float32(rect2)
		dtmp=ntf1d(din.flatten(order='F'),rect1.flatten(order='F'),rect2.flatten(order='F'),niter,n1,verb,rect,ifb,inv,dt,alpha)
	else:
		dtmp=tf1d(din.flatten(order='F'),niter,n1,verb,rect,ifb,inv,dt,alpha)
	#dout: imag,real,basis,w0,dw,nw
	
	if inv==0: #forward transform
		if ifb==1: #output basis functions
			nw=np.int32((dtmp.size-3)/n1/4);
			dout=dtmp[0:n1*nw*2]
			basis=dtmp[n1*nw*2:n1*nw*4]
			w0=dtmp[n1*nw*4]/np.pi/2
			dw=dtmp[n1*nw*4+1]/np.pi/2
			nw2=np.int32(dtmp[n1*nw*4+2])
			if nw2 != nw:
				print('nw=',nw,'nw2=',nw2,'dimension discrepancy')
			else:
				print('dimension consistent')
			return dout,basis,w0,dw,nw
		else:
			nw=np.int32((dtmp.size-3)/n1/2);
			dout=dtmp[0:n1*nw*2]
			w0=dtmp[n1*nw*2]/np.pi/2
			dw=dtmp[n1*nw*2+1]/np.pi/2
			nw2=np.int32(dtmp[n1*nw*2+2])
			print(dtmp.size,w0,dw,nw2)
			if nw2 != nw:
				print('nw=',nw,'nw2=',nw2,'dimension discrepancy')
			else:
				print('dimension consistent')
			return dout,w0,dw,nw
	else:
		return dtmp


	

