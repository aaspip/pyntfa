import pywt #pip install pywavelets
import numpy as np
import matplotlib.pyplot as plt
def spectrum_cwt(t, data, wavelet='mexh', widths=200, cmap='RdBu',
                 colorscale=1, contour_dvision=41, freqmax=100,
                 figsize=(12,3), plot=True,ymin=5, ymax=30, vmax=None):
    '''
    cwt: contineous wavelet transforms for seismic
    '''
    '''
    t: in ms
    data: 1d array
    wavelet:
        ['haar', 'db', 'sym', 'coif', 'bior', 'rbio', 'dmey', 'gaus', 'mexh',
        'morl', 'cgau', 'shan', 'fbsp', 'cmor'] short for:
        ['Haar', 'Daubechies', 'Symlets', 'Coiflets', 'Biorthogonal',
        'Reverse biorthogonal','Discrete Meyer(FIR Approximation)','Gaussian',
        'Mexican hat wavelet', 'Morlet wavelet', 'Complex Gaussian wavelets',
        'Shannon wavelets', 'Frequency B-Spline wavelets',
        'Complex Morlet wavelets']
    widths: the number of frequencies
    colorscale: change color scale, the smaller the darker the color
    contour_division: the number of color divided
    freqmax: the maximum of frequency
    plot: True or False
    '''
    t = np.array(t).reshape(len(t),)
    data = np.array(data).reshape(len(data),)
    
    widths = np.arange(1, widths)
    sampling_rate = 1000 * len(t)/(max(t) - min(t))
    cwtmatr, freqs = pywt.cwt(data, widths, wavelet, 1/sampling_rate)
    
    maximum = colorscale * max(abs(cwtmatr.max()), abs(cwtmatr.min()))
    if vmax != None:  maximum = vmax
    
    cwtmatr_c = np.clip(cwtmatr, -maximum, maximum, cwtmatr)
    levels = np.linspace(-maximum, maximum, contour_dvision)
    
    return cwtmatr, freqs

#import np

from pyseistr import genflat
data=genflat();
n1=data.shape[0]
dt = 0.004
t = np.arange(n1)*dt*1000
cwtmatr, freqs =spectrum_cwt(t, data[:,1], wavelet='morl',
                 widths=51, cmap='RdBu', colorscale=1, contour_dvision=41,
                 freqmax=125, figsize=(12,3), plot=True, ymin=0, ymax=120)

import pyntfa as ntfa
#n1=din.shape[0]
dout,w0,dw,nw = ntfa.ntfa1d(data[:,1],dt=0.004,niter=10,rect=3,ifb=0,inv=0)
dout=dout.reshape([n1,nw,2],order='F');
dtf=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1];
#print(dout[:,:,0].max(),dout[:,:,0].min()) ##should be [0.0864,-0.0876]
## NOTE: dout[:,:,0] is imaginary and dout[:,:,1] is real

## plot TF spectrum
plt.figure();
plt.imshow(dtf,cmap=plt.cm.jet, interpolation='none', extent=[0,nw*dw-dw,n1*1-1,0],aspect='auto');plt.xlabel('Frequency (Hz)');plt.ylabel('Time (s)');plt.title('Time-frequency Spectrum')
plt.show();

t = np.arange(n1)*dt

xmin=0
xmax=120
contour_dvision=41
fig = plt.figure(figsize=(16, 8))
ax=plt.subplot(1,2,1)
maximum = max(abs(cwtmatr.max()), abs(cwtmatr.min()))
#cwtmatr_c = np.clip(cwtmatr, -maximum, maximum, cwtmatr)
levels = np.linspace(-maximum, maximum, contour_dvision)
plt.contourf(freqs, t, cwtmatr.transpose(), levels=51, cmap='RdBu')
#plt.ylim(min(freqs), freqmax) #
#plt.yscale('log')
plt.xlim(xmin, xmax)
#plt.xlim(t.min(), t.max()+1)
plt.ylabel(u"Time (s)")
plt.xlabel(u"Frequency (Hz)")
plt.title(u"Time-frequency Spectrum (CWT)")
plt.colorbar()
plt.tight_layout()


ax=plt.subplot(1,2,2)
maximum = max(abs(dtf.max()), abs(dtf.min()))
levels = np.linspace(-maximum, maximum, contour_dvision)
freqs=np.linspace(0,nw-1,nw)*dw;
plt.contourf(freqs, t, dtf, levels=51, cmap='RdBu')
#plt.ylim(min(freqs), freqmax) #
#plt.yscale('log')
plt.xlim(xmin, xmax)
#plt.xlim(t.min(), t.max()+1)
plt.ylabel(u"Time (s)")
plt.xlabel(u"Frequency(Hz)")
plt.title(u"Time-frequency Spectrum (NTFA)")
plt.colorbar()
plt.tight_layout()
plt.savefig('test_pyntfa_syn1d_noise_cwtVSntfa.png')
plt.show()



## plot noise

from pyseistr import genflat
data=genflat();
trace=data[:,0]
trace=trace/np.max(np.max(trace));
n1=data.shape[0]
dt = 0.004
t = np.arange(n1)*dt
ttmp = np.arange(n1)*dt*1000

fig = plt.figure(figsize=(24, 16))

i=1;
# for sigma in [0,0.2,0.4,0.6,0.8]:
#     np.random.seed(202324+i);
for sigma in [0,0.1,0.2,0.3,0.4]:
    np.random.seed(2023);
    ax=plt.subplot(2,10,2*i-1)
    noise=(np.random.rand(n1)*2-1)*sigma;
    tracen=trace+noise
    plt.plot(tracen,t);plt.xlabel('Amplitude',size=24);
    plt.gca().xaxis.set_tick_params(labelsize=20)
    if i==1:
        plt.ylabel('Time (s)',size=24);plt.gca().yaxis.set_tick_params(labelsize=20)
    else:
        plt.gca().set_yticks([]);

    plt.gca().invert_yaxis();
    
    plt.gca().text(-0.15,1.01,['(a)','(b)','(c)','(d)','(e)'][i-1],transform=plt.gca().transAxes,size=32,weight='normal')
    
    ax=plt.subplot(2,10,2*i)
    cwtmatr, freqs =spectrum_cwt(ttmp, tracen, wavelet='morl',
                 widths=51, cmap='RdBu', colorscale=1, contour_dvision=41,
                 freqmax=125, figsize=(12,3), plot=True, ymin=0, ymax=120)

    maximum = max(abs(cwtmatr.max()), abs(cwtmatr.min()))
    cwtmatr_c = np.clip(cwtmatr, -maximum, maximum, cwtmatr)
    levels = np.linspace(-maximum, maximum, contour_dvision)
    plt.contourf(freqs, t, cwtmatr_c.transpose(), levels=51, cmap='RdBu')
    plt.xlim(xmin, xmax)
    plt.xlabel(u"Frequency (Hz)",size=24)
    plt.gca().xaxis.set_tick_params(labelsize=20)
    plt.tight_layout()
    plt.title('$\sigma$=%g'%sigma,size=24)
    plt.gca().set_yticks([]);
    plt.gca().invert_yaxis();
    
    i=i+1;
# plt.savefig('test_pyntfa_syn1d_noise_cwtnoise.png')
# plt.show()


## Below is NTFA
# fig = plt.figure(figsize=(30, 20))
i=1;
# for sigma in [0,0.2,0.4,0.6,0.8]:
for sigma in [0,0.1,0.2,0.3,0.4]:
    np.random.seed(2023);
#     ax=plt.subplot(2,10,2*i-1+10)
    noise=(np.random.rand(n1)*2-1)*sigma;
    tracen=trace+noise
#     plt.plot(tracen,t);plt.xlabel('Amplitude');
#     if i==1:
#         plt.ylabel('Time (s)');
#     else:
#         plt.gca().set_yticks([]);
        
#     plt.gca().invert_yaxis();
    
#     plt.gca().text(-0.15,1.01,['(f)','(g)','(h)','(i)','(j)'][i-1],transform=plt.gca().transAxes,size=16,weight='normal')
    
    ax=plt.subplot(2,10,2*i+10)
    dout,w0,dw,nw=ntfa.ntfa1d(tracen,dt=0.004,niter=10,rect=3,ifb=0,inv=0)
    dout=dout.reshape([n1,nw,2],order='F');
    dtf=dout[:,:,0]*dout[:,:,0]+dout[:,:,1]*dout[:,:,1];
    freqs=np.linspace(0,nw-1,nw)*dw;
    plt.contourf(freqs, t, dtf, levels=51, cmap='RdBu')
    maximum = max(abs(dtf.max()), abs(dtf.min()))
    dtf_c = np.clip(dtf, -maximum, maximum, dtf)
    levels = np.linspace(-maximum, maximum, contour_dvision)
    plt.contourf(freqs, t, dtf_c, levels=51, cmap='RdBu')
    plt.xlim(xmin, xmax)
    plt.xlabel(u"Frequency (Hz)",size=24)
    plt.gca().xaxis.set_tick_params(labelsize=20)
    plt.tight_layout()
    plt.title('$\sigma$=%g'%sigma,size=24)
    plt.gca().set_yticks([]);
    plt.gca().invert_yaxis();
    
    i=i+1;
plt.savefig('test_pyntfa_noise_cwtVSntfa.png')
plt.savefig('test_pyntfa_noise_cwtVSntfa.pdf')
plt.show()


