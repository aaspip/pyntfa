# Demo of PcP data for an earthquake occured in western Pacific subuduction zone (52.4218,159.9224) recorded by the USArray

# Please download the three files from
# https://github.com/aaspip/data/blob/main/PcP_raw.bin
# https://github.com/aaspip/data/blob/main/PcP_processed.bin
# https://github.com/aaspip/data/blob/main/PcP_noise.bin

import matplotlib.pyplot as plt
import pyntfa as ntfa
import numpy as np

from obspy import read, UTCDateTime

# PcP data example
import numpy as np

# Specify the filename
filename = 'PcP_raw.bin'

# Define the shape of the matrix
# Replace these values with the actual dimensions of your matrix
num_rows = 449   
num_cols = 389

# Calculate the total number of elements
total_elements = num_rows * num_cols

# Read the binary data
with open(filename, 'rb') as file:
    # Read the data as double-precision floats
    data = np.fromfile(file, dtype=np.float64, count=total_elements)
# Reshape the data to the original matrix shape
d0 = data.reshape((num_rows, num_cols), order='F')  # 'F' means column-major

# Specify the filename
filename = 'PcP_processed.bin'
# Read the binary data
with open(filename, 'rb') as file:
    # Read the data as double-precision floats
    data = np.fromfile(file, dtype=np.float64, count=total_elements)
# Reshape the data to the original matrix shape
d1 = data.reshape((num_rows, num_cols), order='F')  # 'F' means column-major

# Specify the filename
filename = 'PcP_noise.bin'
# Read the binary data
with open(filename, 'rb') as file:
    # Read the data as double-precision floats
    data = np.fromfile(file, dtype=np.float64, count=total_elements)
# Reshape the data to the original matrix shape
d2 = data.reshape((num_rows, num_cols), order='F')  # 'F' means column-major

dt = 0.2
earth0 = d0[:,100]
earth1 = d1[:,100]
earth2 = d2[:,100]
# set time indices for different phases
tp = 10
tpcp =  50.3097
times=np.array([tp,tpcp])
phases=['P','PcP']

#plt.plot(earth);
#plt.show()
n1=earth0.shape[0]

#dout,w0,dw,nw=ntfa.ntfa1d(earth,dt=1,niter=10,rect=20,ifb=0,inv=0)
dout0,basis0,w0,dw,nw = ntfa.ntfa1d(earth0,dt=dt,niter=10,rect=10,ifb=1,inv=0)
dout0=dout0.reshape([n1,nw,2],order='F');
basis0=basis0.reshape([n1,nw,2],order='F')

dout1,basis1,w0,dw,nw = ntfa.ntfa1d(earth1,dt=dt,niter=10,rect=10,ifb=1,inv=0)
dout1=dout1.reshape([n1,nw,2],order='F');
basis1=basis1.reshape([n1,nw,2],order='F')

dout2,basis2,w0,dw,nw = ntfa.ntfa1d(earth2,dt=dt,niter=10,rect=10,ifb=1,inv=0)
dout2=dout2.reshape([n1,nw,2],order='F');
basis2=basis2.reshape([n1,nw,2],order='F')

## Time-frequency spectra
dtf0=dout0[:,:,0]*dout0[:,:,0]+dout0[:,:,1]*dout0[:,:,1];
dtf1=dout1[:,:,0]*dout1[:,:,0]+dout1[:,:,1]*dout1[:,:,1];
dtf2=dout2[:,:,0]*dout2[:,:,0]+dout2[:,:,1]*dout2[:,:,1];

freqs=np.linspace(0,nw-1,nw)*dw;
t=np.linspace(0,(n1-1)*dt,n1)

fig = plt.figure(figsize=(12, 6))
plt.subplot(3,2,1)
plt.plot(t,earth0,'k',linewidth=1);plt.ylim(-1,1.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Raw');
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii])
plt.gca().text(-0.18,1,'(a)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.subplot(3,2,3)
plt.plot(t,earth1,'k',linewidth=1);plt.ylim(-1,1.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('Processed');
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii])
plt.gca().text(-0.18,1,'(b)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.subplot(3,2,5)
plt.plot(t,earth2,'k',linewidth=1);plt.ylim(-1,1.2);plt.ylabel('Amplitude');plt.title('Removed noise');
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii])
plt.xlabel('Time (s)');
plt.gca().text(-0.18,1,'(c)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.subplot(3,2,2)
plt.imshow(dtf0.T,clim=(0, 0.0005),cmap=plt.cm.jet, interpolation='none', extent=[0,(n1*1-1)*dt,0,(nw*dw-dw)],origin='lower',aspect='auto');
plt.ylabel('Frequency (Hz)');plt.title('Time-frequency Spectrum')
plt.ylim([0, 1])
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'w--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii],color='w')
plt.gca().text(-0.18,1,'(d)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.subplot(3,2,4)
plt.imshow(dtf1.T,clim=(0, 0.00005),cmap=plt.cm.jet, interpolation='none', extent=[0,(n1*1-1)*dt,0,(nw*dw-dw)],origin='lower',aspect='auto');
plt.ylabel('Frequency (Hz)');
plt.ylim([0, 1])
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'w--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii],color='w')
plt.gca().text(-0.18,1,'(e)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.subplot(3,2,6)
plt.imshow(dtf2.T,clim=(0, 0.0005),cmap=plt.cm.jet, interpolation='none', extent=[0,(n1*1-1)*dt,0,(nw*dw-dw)],origin='lower',aspect='auto');
plt.ylabel('Frequency (Hz)');plt.xlabel('Time (s)');
plt.ylim([0, 1])
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'w--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii],color='w')
plt.gca().text(-0.18,1,'(f)',transform=plt.gca().transAxes,size=18,weight='normal')
plt.savefig('test_pyntfa_pcp.png',format='png',dpi=300)
plt.savefig('test_pyntfa_pcp.pdf',format='pdf',dpi=300)
plt.show();
