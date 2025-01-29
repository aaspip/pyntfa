import pyntfa as ntfa
import numpy as np
import matplotlib.pyplot as plt

din=np.random.randn(1001)
dout=ntfa.tf(din)

err=np.linalg.norm(din-dout)
print(err)

# Download the data from https://github.com/aaspip/data/blob/main/cchirps.bin
fid=open("cchirps.bin","rb");
din = np.fromfile(fid, dtype = np.float32, count = 512).reshape([512,1],order='F')
# plt.plot(din)
# plt.show();

n1=din.shape[0]
dout,w0,dw,nw = ntfa.ntfa1d(din,dt=1,niter=10,rect=7,ifb=0,inv=0,ifn=0)
dout1,basis,w02,dw2,nw2 = ntfa.ntfa1d(din,dt=1,niter=10,rect=7,ifb=1,inv=0,ifn=0)

dout=dout.reshape([n1,nw,2],order='F');
basis=basis.reshape([n1,nw,2],order='F')

print('First version:',dout[:,:,0].max(),dout[:,:,0].min()) ##should be [0.0864,-0.0876]
## NOTE: dout[:,:,0] is imaginary and dout[:,:,1] is real

## Inverse transform
trace=ntfa.ntfa1d(dout,dt=1,niter=10,rect=7,ifb=0,inv=1);

##### Below is for non-stationary Time-frequency regularization

## create radius
r_ref=7 #reference radius
r_max=14 #reference radius

tmp=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1]
tmp=tmp/tmp.max()
tmp[tmp<0.05]=0
tmp[tmp>=0.05]=1
rectn0=tmp*(r_ref-r_max)+r_max 	#first axis radius 	(time)
rectn1=tmp*(-2)+3				#second axis radius	(frequency)

fig = plt.figure(figsize=(16, 8))
plt.subplot(1,2,1)
plt.imshow(rectn0,aspect='auto');
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
plt.title('First Axis (Time) smoothing');

plt.subplot(1,2,2)
plt.imshow(rectn1,aspect='auto');
plt.colorbar(orientation='horizontal',shrink=0.6,label='Traveltime (s)');
plt.title('Second Axis (Frequency) smoothing');

plt.show()

n1=din.shape[0]
doutn,w0,dw,nw = ntfa.ntfa1d(din,dt=1,niter=10,rect=7,ifb=0,inv=0,ifn=1,rect1=rectn0,rect2=rectn1)
doutn1,basisn,w02,dw2,nw2 = ntfa.ntfa1d(din,dt=1,niter=10,rect=7,ifb=1,inv=0,ifn=1,rect1=rectn0,rect2=rectn1)
doutn=doutn.reshape([n1,nw,2],order='F');
basisn=basisn.reshape([n1,nw,2],order='F')
print('Second version:',doutn[:,:,0].max(),doutn[:,:,0].min()) ##should be [0.08177045 -0.081738934]
## Reconstruct trace
tracen=ntfa.ntfa1d(doutn,dt=1,niter=10,rect=7,ifb=0,inv=1,ifn=1,rect1=rectn0,rect2=rectn1);


## Better visualization
dt=1;
t=np.linspace(0,(n1-1)*dt,n1)
fig = plt.figure(figsize=(16, 8))
plt.subplot(3,3,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-3,3);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.gca().text(-0.15,1,'(a)',transform=plt.gca().transAxes,size=20,weight='normal')
plt.subplot(3,3,4)
plt.plot(t,trace,'k',linewidth=1);plt.ylim(-3,3);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.gca().text(-0.15,1,'(b)',transform=plt.gca().transAxes,size=20,weight='normal')

plt.subplot(3,3,7)
plt.plot(t,din.flatten()-trace,'k',linewidth=1);plt.ylim(-3,3);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');
plt.gca().text(-0.15,1,'(c)',transform=plt.gca().transAxes,size=20,weight='normal')


plt.subplot(1,3,2)
plt.imshow(dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1],clim=(0, 0.001),cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum')
plt.gca().text(-0.15,1,'(d)',transform=plt.gca().transAxes,size=20,weight='normal')

plt.subplot(1,3,3)
plt.imshow(basis[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Basis function (Real)')
plt.gca().text(-0.15,1,'(e)',transform=plt.gca().transAxes,size=20,weight='normal')
plt.colorbar(orientation='horizontal',cax=fig.add_axes([0.35,0.05,0.1,0.015]),shrink=1);

plt.savefig('test_pyntfa_syn1d1.png',format='png',dpi=300)
plt.savefig('test_pyntfa_syn1d1.pdf',format='pdf',dpi=300)
plt.show();

fig = plt.figure(figsize=(16, 8))
plt.subplot(3,3,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-3,3);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.gca().text(-0.15,1,'(a)',transform=plt.gca().transAxes,size=20,weight='normal')

plt.subplot(3,3,4)
plt.plot(t,tracen,'k',linewidth=1);plt.ylim(-3,3);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.gca().text(-0.15,1,'(b)',transform=plt.gca().transAxes,size=20,weight='normal')

plt.subplot(3,3,7)
plt.plot(t,din.flatten()-tracen,'k',linewidth=1);plt.ylim(-3,3);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');
plt.gca().text(-0.15,1,'(c)',transform=plt.gca().transAxes,size=20,weight='normal')

plt.subplot(1,3,2)
plt.imshow(doutn[:,:,0]*doutn[:,:,0]+doutn[:,:,1]*doutn[:,:,1],clim=(0, 0.001),cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum')
plt.gca().text(-0.15,1,'(d)',transform=plt.gca().transAxes,size=20,weight='normal')

plt.subplot(1,3,3)
plt.imshow(basisn[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Basis function (Real)')
plt.gca().text(-0.15,1,'(e)',transform=plt.gca().transAxes,size=20,weight='normal')
plt.colorbar(orientation='horizontal',cax=fig.add_axes([0.35,0.05,0.1,0.015]),shrink=1);

plt.savefig('test_pyntfa_syn1d2.png',format='png',dpi=300)
plt.savefig('test_pyntfa_syn1d2.pdf',format='pdf',dpi=300)
plt.show();


## Below is for benchmarking
# from pylib.io import binread
# 
# cchirps=binread('/Users/chenyk/data/datapath/chenyk.rr/tf/cross/cchirps.rsf@',n1=512,n2=1)
# fig = plt.figure(figsize=(16, 8))
# plt.subplot(3,3,1)
# plt.plot(t,din,'k',linewidth=1);plt.ylim(-3,3);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
# plt.subplot(3,3,4)
# plt.plot(t,cchirps,'k',linewidth=1);plt.ylim(-3,3);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Madagascar');
# plt.subplot(3,3,7)
# plt.plot(t,din.flatten()-cchirps,'k',linewidth=1);plt.ylim(-3,3);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Error');
# plt.show()

# ltft_n_real=binread('/Users/chenyk/data/datapath/chenyk.rr/tf/cross/ltft_n_real.rsf@',n1=512,n2=257)
# ltft_n_imag=binread('/Users/chenyk/data/datapath/chenyk.rr/tf/cross/ltft_n_imag.rsf@',n1=512,n2=257)
# ltft_s_real=binread('/Users/chenyk/data/datapath/chenyk.rr/tf/cross/ltft_s_real.rsf@',n1=512,n2=257)
# ltft_s_imag=binread('/Users/chenyk/data/datapath/chenyk.rr/tf/cross/ltft_s_imag.rsf@',n1=512,n2=257)
# 
# fig = plt.figure(figsize=(16, 8))
# plt.imshow(np.concatenate([doutn[:,:,0],ltft_n_imag,doutn[:,:,0]-ltft_n_imag],axis=1));
# plt.show()
# 
# fig = plt.figure(figsize=(16, 8))
# plt.imshow(np.concatenate([doutn[:,:,1],ltft_n_real,doutn[:,:,1]-ltft_n_real],axis=1));
# plt.show()
# 
# fig = plt.figure(figsize=(16, 8))
# plt.imshow(np.concatenate([dout[:,:,0],ltft_s_imag,dout[:,:,0]-ltft_s_imag],axis=1));
# plt.show()
# 
# fig = plt.figure(figsize=(16, 8))
# plt.imshow(np.concatenate([dout[:,:,1],ltft_s_real,dout[:,:,1]-ltft_s_real],axis=1));
# plt.show()
# 



