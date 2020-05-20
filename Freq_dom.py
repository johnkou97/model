import numpy as np
import math

import matplotlib.pyplot as plt
import h5py
def keys(f):
    return [key for key in f.keys()]
import scipy
from scipy import signal
from scipy.fftpack import fft, fftshift ,ifft,rfft,fftfreq,rfftfreq

from scipy.interpolate import CubicSpline as spline
#h5py.run_tests()
c=2.9979e10
G=6.67408e-8
Msun=1.989e33
Length = G*Msun/c**2
Time = Length/c
Frequency=1/Time

def fre_do(x,y,mass):
    fd=fft(y)
    N=len(y)
    T=x[1]-x[0]
    xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))/mass
    fq=fftfreq(len(y))
    mask=fq>=0
    fd=2.0*(fd/N)
    fd=fd[mask]
    fd=abs(fd)

    return xf,fd


def analyze(file,mass):



    rhM=rh[:,1]
    time=rh[:,0]

    peaks,prop=scipy.signal.find_peaks(abs(rhM))
    ampls=rhM[peaks]
    merg=np.amax(abs(ampls))
    merg=np.where(abs(ampls)==merg)
    merg=int(merg[0])
    t0=peaks[merg]

    ampl=rhM[t0:]
    tim=time[t0:]

    tuk=signal.tukey(len(ampl))
    dat=ampl*tuk

    fq,fd=fre_do(tim,dat,mass)

    mx=np.where(fd==np.amax(fd))[0][0]
    freq=fq[mx]
    amp=fd[mx]
    return freq,amp


SLy=1
H4=2
MS1=3
MPA1=4
ALF2=5
MS1b=6
ms1b=6
ENG=7


BHBlp=9
DD2=10
LS220=11
SFHo=12


#read metadata for BAM
mas2=np.zeros(31)
eos2=np.zeros(31)
i=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'metadata/BAM:0%s%s%s.txt' %(m,k,j)
            try:
                f=open(name)
                lines=f.readlines()

                exec(lines[8])
                mas2[i]=id_mass
                if i==0:
                    eos2[i]=8 #by hand because 2H can not be given a value
                if i>0:
                    exec(lines[15])
                    eos2[i]=id_eos


                i=i+1
            except OSError:
                pass





#read metadata for THC
mas1=np.zeros(23)
eos1=np.zeros(23)
i=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'metadata/THC:0%s%s%s.txt' %(m,k,j)
            try:
                f=open(name)
                lines=f.readlines()

                exec(lines[8])
                mas1[i]=id_mass
                exec(lines[15])
                eos1[i]=id_eos


                i=i+1
            except OSError:
                pass

#analyze data for BAM

freq2=np.zeros(31)
amp2=np.zeros(31)
i=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'data/BAM:0%s%s%s.h5' %(m,k,j)
            try:
                file=h5py.File(name,'r')
                dat = list(file["/rh_22"])
                rh = np.array(file["/rh_22/%s" %dat[-1]])

                freq2[i],amp2[i]=analyze(file,mas2[i])

                i=i+1
            except OSError:
                pass


#analyze data for THC


freq1=np.zeros(23)
amp1=np.zeros(23)
i=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'data/THC:0%s%s%s.h5' %(m,k,j)
            try:
                file=h5py.File(name,'r')
                dat = list(file["/rh_22"])
                rh = np.array(file["/rh_22/Rh_l2_m2_r00400.txt" ])

                freq1[i],amp1[i]=analyze(file,mas1[i])

                i=i+1
            except OSError:
                pass




freq=np.concatenate((freq1,freq2))
mas=np.concatenate((mas1,mas2))
eos=np.concatenate((eos1,eos2))
amp=np.concatenate((amp1,amp2))

#mtov for different eos

mtov=np.zeros(len(eos))
for i in range(len(eos)):
    if eos[i]==1:
        mtov[i]=2.06
    if eos[i]==2:
        mtov[i]=2.03
    if eos[i]==3:
        mtov[i]=2.77
    if eos[i]==4:
        mtov[i]=2.77
    if eos[i]==5:
        mtov[i]=1.99
    if eos[i]==6:
        mtov[i]=2.76
    if eos[i]==7:
        mtov[i]=2.25
    if eos[i]==8:
        mtov[i]=2.83
    if eos[i]==9:
        mtov[i]=2.10
    if eos[i]==10:
        mtov[i]=2.42
    if eos[i]==11:
        mtov[i]=2.04
    if eos[i]==12:
        mtov[i]=2.06


lamda=np.load('results/lamda.npy')
k=(3/16)*lamda
a=-131.7010
zeta=k+a*mas/mtov


np.save('results/mas.npy',mas)
np.save('results/freq.npy',freq)
np.save('results/zeta.npy',zeta)
np.save('results/amp.npy',amp)
