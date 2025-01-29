import matplotlib.pyplot as plt
import pyntfa as ntfa
from pyntfa import binread
import numpy as np

# Download data from https://github.com/aaspip/data/blob/main/usarray_trace.bin
earth=binread('usarray_trace.bin',n1=5400,n2=1)

dt=1

# set time indices for different phases
tp=790
tpp=1050
ts=1550
tss=2100
tr=2600

times=np.array([tp,tpp,ts,tss,tr])+586.74
phases=['P','PP','S','SS','Rayleigh']

#plt.plot(earth);
#plt.show()
n1=5400;
#dout,w0,dw,nw=ntfa.ntfa1d(earth,dt=1,niter=10,rect=20,ifb=0,inv=0)
dout,basis,w0,dw,nw = ntfa.ntfa1d(earth,dt=1,niter=10,rect=20,ifb=1,inv=0)
dout=dout.reshape([n1,nw,2],order='F');
basis=basis.reshape([n1,nw,2],order='F')

## Reconstruct trace
earth2=ntfa.ntfa1d(dout,dt=1,niter=10,rect=7,ifb=0,inv=1);

## Time-frequency spectra
dtf=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1];
freqs=np.linspace(0,nw-1,nw)*dw;

dt=1;
t=np.linspace(0,(n1-1)*dt,n1)
fig = plt.figure(figsize=(16, 8))
plt.subplot(3,3,1)
plt.plot(t,earth,'k',linewidth=1);plt.ylim(-0.2,0.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii],ymax-(ymax-ymin)*0.2,phases[ii])

plt.subplot(3,3,4)
plt.plot(t,earth2,'k',linewidth=1);plt.ylim(-0.2,0.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii],ymax-(ymax-ymin)*0.2,phases[ii])
    
plt.subplot(3,3,7)
plt.plot(t,earth.flatten()-earth2,'k',linewidth=1);plt.ylim(-0.2,0.2);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');

plt.subplot(1,3,2)
plt.imshow(dtf,clim=(0, 0.00000001),cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum')
#adding labels
xmin, xmax = plt.gca().get_xlim()
for ii in range(len(times)):
    plt.plot([xmin,xmax],[times[ii],times[ii]],'w--',linewidth=1);
    plt.text(xmax-(xmax-xmin)*0.2,times[ii],phases[ii],color='w')

plt.subplot(1,3,3)
plt.imshow(basis[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Basis function (Real)')

plt.savefig('test_pyntfa_usarray1d1.png',format='png',dpi=300)
plt.show();



## create radius
r_ref=20 #reference radius
r_max=40 #reference radius

tmp=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1]
tmp=tmp/tmp.max()
tmp[tmp<0.002]=0
tmp[tmp>=0.002]=1
rectn0=tmp*(r_ref-r_max)+r_max     #first axis radius     (time)
rectn1=tmp*(-2)+3                #second axis radius    (frequency)

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

n1=earth.shape[0]
#doutn,w0,dw,nw = ntfa.ntfa1d(earth,dt=1,niter=10,rect=7,ifb=0,inv=0,ifn=1,rect1=rectn0,rect2=rectn1)
doutn,basisn,w02,dw2,nw2 = ntfa.ntfa1d(earth,dt=1,niter=10,rect=7,ifb=1,inv=0,ifn=1,rect1=rectn0,rect2=rectn1)
doutn=doutn.reshape([n1,nw,2],order='F');
basisn=basisn.reshape([n1,nw,2],order='F')
print('Second version:',doutn[:,:,0].max(),doutn[:,:,0].min()) ##should be [0.08177045 -0.081738934]
## Reconstruct trace
earthn=ntfa.ntfa1d(doutn,dt=1,niter=10,rect=7,ifb=0,inv=1,ifn=1,rect1=rectn0,rect2=rectn1);
dtfn=doutn[:,:,0]*doutn[:,:,0]+doutn[:,:,1]*doutn[:,:,1]

## Better visualization
dt=1;
t=np.linspace(0,(n1-1)*dt,n1)
fig = plt.figure(figsize=(16, 8))
plt.subplot(3,3,1)
plt.plot(t,earth,'k',linewidth=1);plt.ylim(-0.2,0.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Input');
plt.gca().text(-0.18,1,'(a)',transform=plt.gca().transAxes,size=18,weight='normal')

#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii],ymax-(ymax-ymin)*0.2,phases[ii])
    
plt.subplot(3,3,4)
plt.plot(t,earthn,'k',linewidth=1);plt.ylim(-0.2,0.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Reconstruction');
plt.gca().text(-0.18,1,'(b)',transform=plt.gca().transAxes,size=18,weight='normal')

#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii],ymax-(ymax-ymin)*0.2,phases[ii])
    
plt.subplot(3,3,7)
plt.plot(t,earth.flatten()-earthn,'k',linewidth=1);plt.ylim(-0.2,0.2);plt.ylabel('Amplitude');plt.xlabel('Time (s)');plt.title('Reconstruction Error');
plt.gca().text(-0.18,1,'(c)',transform=plt.gca().transAxes,size=18,weight='normal')


plt.subplot(1,3,2)
plt.imshow(dtfn,clim=(0, 0.00000001),cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum')
plt.gca().text(-0.18,1,'(d)',transform=plt.gca().transAxes,size=18,weight='normal')
#adding labels
xmin, xmax = plt.gca().get_xlim()
for ii in range(len(times)):
    plt.plot([xmin,xmax],[times[ii],times[ii]],'w--',linewidth=1);
    plt.text(xmax-(xmax-xmin)*0.2,times[ii],phases[ii],color='w')

plt.subplot(1,3,3)
plt.imshow(basis[:,:,1],cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Basis function (Real)')
plt.gca().text(-0.18,1,'(e)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.colorbar(orientation='horizontal',cax=fig.add_axes([0.35,0.05,0.1,0.015]),shrink=1);

plt.savefig('test_pyntfa_usarray1d.png',format='png',dpi=300)
plt.show();




