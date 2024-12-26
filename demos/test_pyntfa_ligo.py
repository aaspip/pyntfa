import numpy as np
import matplotlib.pyplot as plt
from pyntfa import ntfa1d
from pynpre import stft1d,st1d


# Data loading
# Download the data from https://github.com/aaspip/data/blob/main/ligosig.bin
fid=open("ligosig.bin","rb");
din = np.fromfile(fid, dtype = np.float32, count = 2260).reshape([2260,1],order='F')
plt.plot(din)
plt.show();


n1=din.shape[0]
dt=0.000244141;
ot=-0.308695

## STFT
dout,w0,dw,nw=stft1d(din,dt=dt,inv=0,opt=1,sym=0,ntw=101,ot=0,wind=0,verb=1)
dout=dout.reshape([n1,nw,2],order='F');
print("dout.shape is ",dout.shape)
## Inverse transform
trace=stft1d(dout,dt=0.000244141,inv=1,opt=1,sym=0,ntw=101,ot=0,wind=0,verb=1)

## Visualization of STFT
t=np.linspace(ot,(n1-1)*dt+ot,n1)
fig = plt.figure(figsize=(12, 8))
plt.subplot(3,2,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.subplot(3,2,3)
plt.plot(t,trace,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.subplot(3,2,5)
plt.plot(t,din.flatten()-trace,'k',linewidth=1);plt.ylim(-4,4);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');

plt.subplot(1,2,2)#clim=(0, 10)
plt.imshow(dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum (STFT)')
plt.savefig('test_pyntfa_stltft_stft_ligo.png',format='png',dpi=300)
plt.show();


## ST
dout,w0,dw,nw=st1d(din,dt=dt,inv=0,flo=0,fhi=0.5,verb=1)
dout=dout.reshape([n1,nw,2],order='F');
print("doutkk.shape is ",dout.shape)
## Inverse transform
trace=st1d(dout,dt=dt,inv=1,flo=0,fhi=0.5,verb=1)

## Visualization of ST
t=np.linspace(ot,(n1-1)*dt+ot,n1)
fig = plt.figure(figsize=(12, 8))
plt.subplot(3,2,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.subplot(3,2,3)
plt.plot(t,trace,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.subplot(3,2,5)
plt.plot(t,din.flatten()-trace,'k',linewidth=1);plt.ylim(-4,4);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');

plt.subplot(1,2,2)#clim=(0, 10)
plt.imshow(dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum (ST)')
plt.savefig('test_pyntfa_stltft_st_ligo.png',format='png',dpi=300)
plt.show();


## NTFA
##### Below is for non-stationary data regularizaton
dout,w0,dw,nw = ntfa1d(din,dt=dt,niter=10,rect=25,ifb=0,inv=0,ifn=0)
dout=dout.reshape([n1,nw,2],order='F');

print('First version:',dout[:,:,0].max(),dout[:,:,0].min()) 
## NOTE: dout[:,:,0] is imaginary and dout[:,:,1] is real

## Visualization of NTFA
t=np.linspace(ot,(n1-1)*dt+ot,n1)
fig = plt.figure(figsize=(12, 8))
plt.subplot(3,2,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.subplot(3,2,3)
plt.plot(t,trace,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.subplot(3,2,5)
plt.plot(t,din.flatten()-trace,'k',linewidth=1);plt.ylim(-4,4);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');

plt.subplot(1,2,2)#clim=(0, 10)
plt.imshow(dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum (NTFA)')
plt.savefig('test_pyntfa_ligo_1.png',format='png',dpi=300)
plt.show();


##### Below is for non-stationary data+model regularizaton
# create radius
r_ref=25 #reference radius
r_max=40 #reference radius

tmp=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1]
tmp=tmp/tmp.max()
tmp[tmp<0.1]=0
tmp[tmp>=0.1]=1
rectn0=tmp*(r_ref-r_max)+r_max     #first axis radius     (time)
rectn1=tmp*(-2)+1                #second axis radius    (frequency)

fig = plt.figure(figsize=(16, 8))
plt.subplot(1,2,1)
plt.imshow(rectn0,aspect='auto');
plt.colorbar(orientation='horizontal',shrink=0.6,label='Radius (samples)');
plt.title('First Axis (Time) smoothing');

plt.subplot(1,2,2)
plt.imshow(rectn1,aspect='auto');
plt.colorbar(orientation='horizontal',shrink=0.6,label='Radius (samples)');
plt.title('Second Axis (Frequency) smoothing');
plt.show()

n1=din.shape[0]
doutn,w0,dw,nw = ntfa1d(din,dt=dt,niter=10,rect=25,ifb=0,inv=0,ifn=1,verb=0,rect1=rectn0,rect2=rectn1)
doutn=doutn.reshape([n1,nw,2],order='F');
print('Second version:',doutn[:,:,0].max(),doutn[:,:,0].min()) 
## Reconstruct trace
tracen=ntfa1d(doutn,dt=0.00024414,niter=10,rect=25,ifb=0,inv=1,ifn=1,rect1=rectn0,rect2=rectn1);

## Visualization of NTFA2
t=np.linspace(ot,(n1-1)*dt+ot,n1)
fig = plt.figure(figsize=(12, 8))
plt.subplot(3,2,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.subplot(3,2,3)
plt.plot(t,tracen,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.subplot(3,2,5)
plt.plot(t,din.flatten()-tracen,'k',linewidth=1);plt.ylim(-4,4);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');

plt.subplot(1,2,2)#clim=(0, 10)
plt.imshow(doutn[:,:,0]*doutn[:,:,0]+doutn[:,:,1]*doutn[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum (NTFA)')
plt.savefig('test_pyntfa_ligo_2.png',format='png',dpi=300)
plt.show();


fig = plt.figure(figsize=(16, 6.4))
plt.subplot(3,3,1)
plt.plot(t,din,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.gca().text(-0.15,1,'(a)',transform=plt.gca().transAxes,size=16,weight='normal')

plt.subplot(3,3,4)
plt.plot(t,tracen,'k',linewidth=1);plt.ylim(-4,4);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.gca().text(-0.15,1,'(b)',transform=plt.gca().transAxes,size=16,weight='normal')

plt.subplot(3,3,7)
plt.plot(t,din.flatten()-tracen,'k',linewidth=1);plt.ylim(-4,4);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');
plt.gca().text(-0.15,1,'(c)',transform=plt.gca().transAxes,size=16,weight='normal')

plt.subplot(1,5,3)#clim=(0, 10)
plt.imshow(dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)',fontsize=10);
plt.title('Time-frequency Spectrum (NTFA)',fontsize=10)
plt.gca().set_xlim(xmin=0,xmax=500);
plt.gca().text(-0.15,1.05,'(d)',transform=plt.gca().transAxes,size=16,weight='normal')

plt.subplot(1,5,4)#clim=(0, 10)
doutstft,w0,dw,nw=stft1d(din,dt=dt,inv=0,opt=1,sym=0,ntw=101,ot=0,wind=0,verb=1)
doutstft=doutstft.reshape([n1,nw,2],order='F');
plt.imshow(doutstft[:,:,0]*doutstft[:,:,0]+doutstft[:,:,1]*doutstft[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');#plt.ylabel('Time (s)');
plt.title('Time-frequency Spectrum (STFT)',fontsize=10)
plt.gca().set_xlim(xmin=0,xmax=500);
plt.gca().text(-0.15,1.05,'(e)',transform=plt.gca().transAxes,size=16,weight='normal')

plt.subplot(1,5,5)#clim=(0, 10)
doutS,w0,dw,nw=st1d(din,dt=dt,inv=0,flo=0,fhi=0.5,verb=1)
doutS=doutS.reshape([n1,nw,2],order='F');
plt.imshow(doutS[:,:,0]*doutS[:,:,0]+doutS[:,:,1]*doutS[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,(n1-1)*dt+ot,ot],aspect='auto');plt.xlabel('Frequency (Hz)');#plt.ylabel('Time (s)');
plt.title('Time-frequency Spectrum (ST)',fontsize=10)
plt.gca().set_xlim(xmin=0,xmax=500);
plt.gca().text(-0.15,1.05,'(f)',transform=plt.gca().transAxes,size=16,weight='normal')

plt.savefig('test_pyntfa_ligo.png',format='png',dpi=300)
plt.show();


