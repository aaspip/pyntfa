# Demo of SS precursors for an Magnitude 7.6 earthquake occured near Tonga (-22.1220, 175.1630) recorded by Station OBN in Russia (55.1146, 36.5674)
import matplotlib.pyplot as plt
import pyntfa as ntfa
import numpy as np

from obspy import read, UTCDateTime

# PcP data example
import numpy as np

# Please download the three files from
# https://github.com/aaspip/data/blob/main/ntfa_ss.bin

# Specify the filename
filename = 'ntfa_ss.bin'

# Define the shape of the matrix
# Replace these values with the actual dimensions of your matrix
num_rows = 1301   
num_cols = 2

# Calculate the total number of elements
total_elements = num_rows * num_cols

# Read the binary data
with open(filename, 'rb') as file:
    # Read the data as double-precision floats
    data = np.fromfile(file, dtype=np.float64, count=total_elements)
# Reshape the data to the original matrix shape
d = data.reshape((num_rows, num_cols), order='F')[:,1]  # 'F' means column-major
t = data.reshape((num_rows, num_cols), order='F')[:,0]  # 'F' means column-major
d = d/np.max(d)

# set time indices for different phases
tss = 0
t410 = -155
t660 = -236
times=np.array([0,t410,t660])
phases=['SS','S410S','S660S']

n1=d.shape[0]

dout,basis,w0,dw,nw = ntfa.ntfa1d(d,dt=1,niter=10,rect=6,ifb=1,inv=0)
dout=dout.reshape([n1,nw,2],order='F');
basis=basis.reshape([n1,nw,2],order='F')

## Time-frequency spectra
dtf=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1];

freqs=np.linspace(0,nw-1,nw)*dw;

fig = plt.figure(figsize=(7, 8))
plt.subplot(4,1,1)
plt.plot(t,d,'k',linewidth=1);
plt.ylim(-1,1.2);plt.gca().set_xticks([]);plt.ylabel('Amplitude');plt.title('SS precursors');
plt.xlim(-400,200)
#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'k--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii])
plt.gca().text(-0.15,1,'(a)',transform=plt.gca().transAxes,size=18,weight='normal')

plt.subplot(4,1,2)
plt.imshow(np.log(dtf.T),clim=(-25,-10),cmap=plt.cm.jet, interpolation='none', extent=[t[0],t[-1],0,(nw*dw-dw)],origin='lower',aspect='auto');
plt.ylabel('Frequency (Hz)');plt.title('Time-frequency Spectrum (NTFA)')
plt.xlim(-400,200)
plt.ylim([0, .2])
# plt.xlabel('Time (s)');
plt.gca().set_xticks([]);

#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'w--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii],color='w')
plt.gca().text(-0.15,1,'(b)',transform=plt.gca().transAxes,size=18,weight='normal')


## For STFT
from pynpre import stft1d,st1d
dout0,w0,dw,nw = stft1d(d,dt=1,inv=0,opt=1,sym=0,ntw=51,ot=0,wind=0,verb=1)
dout=dout0.reshape([n1,nw,2],order='F');
dtf=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1];
plt.subplot(4,1,3)
plt.imshow(np.log(dtf.T),clim=(-2,3),cmap=plt.cm.jet, interpolation='none', extent=[t[0],t[-1],0,(nw*dw-dw)],origin='lower',aspect='auto');
plt.ylabel('Frequency (Hz)');plt.title('Time-frequency Spectrum (STFT)')
plt.xlim(-400,200)
plt.ylim([0, .2])
# plt.xlabel('Time (s)');
plt.gca().set_xticks([]);

#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'w--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii],color='w')
plt.gca().text(-0.15,1,'(c)',transform=plt.gca().transAxes,size=18,weight='normal')


## For ST
dout0,w0,dw,nw = st1d(d,dt=1,inv=0,flo=0,fhi=0.5,verb=1)
dout=dout0.reshape([n1,nw,2],order='F');
dtf=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1];
plt.subplot(4,1,4)
plt.imshow(np.log(dtf.T),clim=(-10,-2),cmap=plt.cm.jet, interpolation='none', extent=[t[0],t[-1],0,(nw*dw-dw)],origin='lower',aspect='auto');
plt.ylabel('Frequency (Hz)');plt.title('Time-frequency Spectrum (ST)')
plt.xlim(-400,200)
plt.ylim([0, .2])
plt.xlabel('Time (s)');

#adding labels
ymin, ymax = plt.gca().get_ylim()
for ii in range(len(times)):
    plt.plot([times[ii],times[ii]],[ymin,ymax],'w--',linewidth=1);
    plt.text(times[ii]+1,ymax-(ymax-ymin)*0.2,phases[ii],color='w')
plt.gca().text(-0.15,1,'(d)',transform=plt.gca().transAxes,size=18,weight='normal')


plt.savefig('test_pyntfa_ss.png',format='png',dpi=300)
plt.savefig('test_pyntfa_ss.pdf',format='pdf',dpi=300)
plt.show();




