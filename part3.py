# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 18:38:32 2018

@author: Dixit Khosla
"""
from scipy.special import expn
import matplotlib.pyplot as plt
import numpy as np
import soundfile 
import math
#import signal 
from scipy import signal


def logmmse(x, Srate, noise_frames=100, Slen=0, eta=4.5, saved_params=None):
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
   
    ksi_min = 10 ** (-25 / 10)

    for k in range(0, Nframes*len2, len2):
        insign = win * x[k:k + Slen]

        spec = np.fft.fft(insign, nFFT, axis=0)
        sig = np.absolute(spec)
        sig2 = sig ** 2

        gammak = np.minimum(sig2 / noise_mu2, 40)

        if Xk_prev.all() == 0:
            #print('love')
            ksi = aa + (1 - aa) * np.maximum(gammak - 1, 0)
        else:
            #print('hate')
            ksi = aa * Xk_prev / noise_mu2 + (1 - aa) * np.maximum(gammak - 1, 0)
            ksi = np.maximum(ksi_min, ksi)

        log_sigma_k = gammak * ksi/(1 + ksi) - np.log(1 + ksi)
        vad_decision = np.sum(log_sigma_k)/Slen
        #print(vad_decision)
        if (vad_decision < eta):
            noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2
           # print('love')

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






x, fs1 = soundfile.read('AfterBeamforming.wav')

output, q= logmmse(x,fs1,noise_frames=100,Slen=320,eta=30)#,saved_params={'noise_mu2': 100)#,'Xk_prev':None,'x_old':None})
soundfile.write('final1.wav',output,fs1) 
#soundfile.write('noise.wav',q['noise_mu2'],fs1) 
print ('Done')