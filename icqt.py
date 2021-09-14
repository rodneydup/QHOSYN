from librosa.core.constantq import icqt
from librosa.core.constantq import cqt
from scipy.io import wavfile
import csv
import pandas as pd
import numpy as np


with open('pre-ifft.csv', 'r') as inputFile:
    # sr, audioIn = wavfile.read('xin.wav')
    # audioFloat = audioIn / 32768.0
    # print(np.shape(audioFloat))
    # analyze = cqt(audioFloat[:,0],sr,n_bins=512,fmin=10,bins_per_octave=64,hop_length=256)
    # print(np.shape(analyze))

    complexData = np.genfromtxt(inputFile, dtype=complex)
    
    complexData = np.reshape(complexData,(256,-1))
    print(np.shape(complexData))
    result = icqt(complexData,sr=48000,hop_length=128,fmin=30,bins_per_octave=64,sparsity=0.99,norm=None,scale=False)
    wavfile.write('out.wav',48000,result)

