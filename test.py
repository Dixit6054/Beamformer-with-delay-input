clear
import numpy as np
import matplotlib.pyplot as plt
from scipy import io
# from scipy import signal
from numpy import correlate
import scipy.io.wavfile
# from cmath import exp
# from scipy import fftpack
rate, sig1 = io.wavfile.read('mixture1.wav')
rate, sig2 = io.wavfile.read('mixture2.wav')
rate, sig3 = io.wavfile.read('mixture3.wav')


# sig1 = signal.decimate(sig1, 4)
# sig1 = signal.decimate(sig1, 4)
# sig2 = signal.decimate(sig2, 4)
# sig2 = signal.decimate(sig2, 4)
# sig3 = signal.decimate(sig3, 4)
# sig3 = signal.decimate(sig3, 4)

crossCor12 = correlate(sig1, sig2, 'full')
# print(np.argmax(crossCor13), len(crossCor13))
# print(len(sig1))
# print(rate)
# crossCor13 = np.correlate(sig1, sig3, "full")
# crossCor23 = np.correlate(sig2, sig3, "full")

# sig1F = fftpack.fft(sig1)
# sig2F = fftpack.fft(sig2)
# sig3F = fftpack.fft(sig3)
#
# print(len(sig1), len(sig2), len(sig3))
#
# crossCorF = sig1F * np.conj(sig2F)
# crossCor = fftpack.ifft(crossCorF)
# temp = np.argsort(crossCor)
# print(temp[-1]/16000, temp[-2]/16000)
# print(crossCor[temp[-1]], crossCor[temp[-2]], crossCor[temp[-3]])
#
# print(len(crossCorF))
# plt.plot(sig1)
# plt.plot(sig2)
# plt.plot(sig3)
# plt.plot(crossCorF)
# plt.show()
# peak12 = np.argsort(crossCor12)
peak12 = np.argsort(-crossCor12)
# # peak23 = np.argsort(crossCor23)
# print(max(crossCor13))
# print(np.argmax(crossCor13))
# print(crossCor13[np.argmax(crossCor13)])
# print(crossCor13[peak13[0]])
newCor = np.zeros(49)
i = 0
while i < len(peak12):
    if 191975 < peak12[i] < 192025:
        print(peak12[i] - 192000, crossCor12[peak12[i]])
        newCor[peak12[i] - 191976] = crossCor12[peak12[i]]
    i += 1

# for i in range(0, np.size(newCor)):
#     if newCor[i] < 0:
#         newCor[i] = 0
# signew1 = sig1
# signew2 = sig2
# signew3 = sig3
#
# for i in range(32000, 191989):
#     signew2[i] = sig2[i+11]
#
# for i in range(32000, 191966):
#     signew3[i] = sig3[i+34]
#
# newSig1 = 1/3 * (sig1 + signew2 + signew3)
#
# signew1 = sig1
# signew2 = sig2
# signew3 = sig3
#
# for i in range(32000, 191989):
#     signew1[i] = sig1[i+11]
#
# for i in range(32000, 191991):
#     signew3[i] = sig3[i+9]
#
# newSig2 = 1/3 * (sig2 + signew1 + signew3)
#
# signew1 = sig1
# signew2 = sig2
# signew3 = sig3
#
# for i in range(32000, 191962):
#     signew1[i] = sig1[i+38]
#
# for i in range(32000, 191982):
#     signew2[i] = sig2[i+18]
#
# newSig3 = 1/3 * (sig3 + signew1 + signew2)
# # # dtheta = np.transpose([exp(-1j*0.0006875), exp(1j*0.00125), exp(1j*0.0006875)])
# # #
# # # sig1F = fftpack.fft(sig1)
# # # sig2F = fftpack.fft(sig2)
# # # sig3F = fftpack.fft(sig3)
# # #
# # # sig = [sig1F, sig3F, sig2F]
# # # newSigF = np.matmul(dtheta, sig)
# # # newSig = fftpack.ifft(newSigF)
# io.wavfile.write("source1.wav", rate, newSig1.astype(sig1.dtype))
# io.wavfile.write("source2.wav", rate, newSig2.astype(sig1.dtype))
# io.wavfile.write("source3.wav", rate, newSig3.astype(sig1.dtype))

# print(np.size(newSig))
#
# plt.figure()
# # # plt.subplot(321)
# plt.plot(newSig)
#
# # plt.subplot(323)
# # plt.plot(crossCor13)
# #
# # plt.subplot(325)
# # plt.plot(crossCor23)
# #
# # plt.subplot(322)
# # plt.plot(crossCor12[11500:12500])
# #
# # plt.subplot(324)
# # plt.plot(crossCor13[11500:12500])
# #
# # plt.subplot(326)
# # plt.plot(crossCor23[11500:12500])
#
# plt.show()
