# -*- coding: utf-8 -*-
"""
Created on Sat Apr  7 21:54:55 2018
@author: Dixit Khosla
"""
from scipy.special import expn
import matplotlib.pyplot as plt
import numpy as np
import soundfile 
import math
#import signal 
from scipy import signal

def logmmse(x, Srate, noise_frames=6, Slen=0, eta=0, saved_params=None):
    print('I am running')
    if Slen == 0:
        Slen = int(math.floor(0.02 * Srate))

    if Slen % 2 == 1:
        Slen = Slen + 1

    PERC = 50
    len1 = math.floor(Slen * PERC / 100)
    len2 = int(Slen - len1)

    win = np.hanning(Slen)
    win = win * len2 / np.sum(win)
    nFFT = 2 * Slen

    x_old = np.zeros(len1)
    Xk_prev = np.zeros(len1)
    Nframes = int(math.floor(len(x) / len2) - math.floor(Slen / len2))
    xfinal = np.zeros(Nframes * len2)

    if saved_params == None:
        noise_mean = np.zeros(nFFT)
        for j in range(0, Slen*noise_frames, Slen):
            noise_mean = noise_mean + np.absolute(np.fft.fft(win * x[j:j + Slen], nFFT, axis=0))
        noise_mu2 = noise_mean / noise_frames ** 2
    else:
        noise_mu2 = saved_params['noise_mu2']
        Xk_prev = saved_params['Xk_prev']
        x_old = saved_params['x_old']

    aa = 0.98
    mu = 0.98
    #eta = 4.5
    ksi_min = 10 ** (-25 / 10)

    for k in range(0, Nframes*len2, len2):
        insign = win * x[k:k + Slen]

        spec = np.fft.fft(insign, nFFT, axis=0)
        sig = np.absolute(spec)
        sig2 = sig ** 2

        gammak = np.minimum(sig2 / noise_mu2, 40)

        if Xk_prev.all() == 0:
            ksi = aa + (1 - aa) * np.maximum(gammak - 1, 0)
        else:
            ksi = aa * Xk_prev / noise_mu2 + (1 - aa) * np.maximum(gammak - 1, 0)
            ksi = np.maximum(ksi_min, ksi)

        log_sigma_k = gammak * ksi/(1 + ksi) - np.log(1 + ksi)
        vad_decision = np.sum(log_sigma_k)/Slen
        if (vad_decision < eta):
            noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2

        A = ksi / (1 + ksi)
        vk = A * gammak
        ei_vk = 0.5 * expn(1, vk)
        hw = A * np.exp(ei_vk)

        sig = sig * hw
        Xk_prev = sig ** 2
        xi_w = np.fft.ifft(hw * spec, nFFT, axis=0)
        xi_w = np.real(xi_w)

        xfinal[k:k + len2] = x_old + xi_w[0:len1]
        x_old = xi_w[len1:Slen]

    return xfinal, {'noise_mu2': noise_mu2, 'Xk_prev': Xk_prev, 'x_old': x_old}

# reading sound files
# sig1, sig2 and sig3 contains data from audio files
# fs1, fs2 and fs3: sampling frequency
# N = length of the audio signal
sig1, fs1 = soundfile.read('mixture1.wav')
sig2, fs2 = soundfile.read('mixture2.wav')
sig3, fs3 = soundfile.read('mixture3.wav')
N = sig1.shape[0]
#print(N)
C13 = np.correlate(sig3,sig1,'full')

#We will plot the C13 for -50 to 50 samples only as the maximum lag can only be 1.5ms
# which is equal to 24 samples according to the given sampling frequency
t = [(value/fs1)*1000.0 for value in range(-50,50)]  # time axis in ms
C13 = [C13[(N-1)+t1] for t1 in range(-50,50)]        # amplitude axis

#cross correlation plot
plt.figure(4)
plt.title('Cross-correlation over -50 to 50 samples ')
plt.plot(t,C13)
plt.xlabel('Time delay (ms)-------------------->')
plt.ylabel('Cross-Correlation of 1st and 3rd microphone----------------->')
plt.show()


t1 = (-1.4375*0.001)/2
t2 = (0.25*0.001)/2
t3 = (1.4375*0.001)/2
N =2048
fs=16000

m1_t = np.mat(np.zeros((3,int((N/2)+1)),dtype=complex))
m2_t = np.mat(np.zeros((3,int((N/2)+1)),dtype=complex))
m3_t = np.mat(np.zeros((3,int((N/2)+1)),dtype=complex))

#steering vectors samples [0,N/2]---> stored in temp variables
for i in range(0,(int(N/2)+1)):
    m1_t[:,i] = (np.mat([1,np.exp((-1j*2*math.pi*fs*i*t1)/N),np.exp((-1j*2*math.pi*fs*i*2*t1)/N)])).T    
    m2_t[:,i] = (np.mat([1,np.exp((-1j*2*math.pi*fs*i*t2)/N),np.exp((-1j*2*math.pi*fs*i*2*t2)/N)])).T
    m3_t[:,i] = (np.mat([1,np.exp((-1j*2*math.pi*fs*i*t3)/N),np.exp((-1j*2*math.pi*fs*i*2*t3)/N)])).T

m1 = np.mat(np.zeros((3,N),dtype=complex))
m2 = np.mat(np.zeros((3,N),dtype=complex))
m3 = np.mat(np.zeros((3,N),dtype=complex))
W = np.mat(np.zeros((3,N),dtype=complex))

for i in range(0,(int(N/2)+1)):
    m1[:,i] = m1_t[:,i]
    m2[:,i] = m2_t[:,i]
    m3[:,i] = m3_t[:,i]

# steering vectors--->[1, N/2 -1]--->filp and conjugate---> placed at first[N/2 + 1,N-1]
for i in range((int(N/2)+1),int(N)):    
    m1[:,i] = np.conj(m1_t[:,int(N-i)])
    m2[:,i] = np.conj(m2_t[:,int(N-i)])
    m3[:,i] = np.conj(m3_t[:,int(N-i)])
alpha = 0.9
for i in range(0,int(N)):
    R = (m1[:,i]*m1[:,i].H) + (m3[:,i]*m3[:,i].H) + ((1-alpha)*np.eye(3))
    W[:,i] = ((R.I)*m2[:,i])/((m2[:,i].H) * (R.I) * (m2[:,i]))
    
W1 =np.conjugate(W[0,:])
W2 =np.conjugate(W[1,:])
W3 =np.conjugate(W[2,:])
# ifft 
h1 =np.fft.ifftshift(np.fft.ifft(W1))  
h2 =np.fft.ifftshift(np.fft.ifft(W2)) 
h3 =np.fft.ifftshift(np.fft.ifft(W3))
y = (signal.lfilter(h1[0,:].real, [1.0], sig1) )+ (signal.lfilter(h2[0,:].real, [1.0], sig2)) +( signal.lfilter(h3[0,:].real, [1.0], sig3))
soundfile.write('AfterBeamforming.wav',y,fs1) 
#part 3, de-noising
output, q= logmmse(y,fs1, noise_frames=100,Slen=320,eta=29.4)
soundfile.write('final.wav',output,fs1) 