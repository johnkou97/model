import numpy as np
import math
import matplotlib.pyplot as plt
import h5py
import os
def keys(f):
    return [key for key in f.keys()]
import scipy
from scipy import signal
from scipy.fftpack import fft, fftshift ,ifft
import gwpy
from gwpy.timeseries import TimeSeries
from scipy.interpolate import CubicSpline as spline
#h5py.run_tests()
c=2.9979e10
G=6.67408e-8
Msun=1.989e33
Length = G*Msun/c**2
Time = Length/c


def analyze(rh):
    peaks,prop=scipy.signal.find_peaks(abs(rh[:,1]))
    mx=np.where(abs(rh[peaks])==np.amax(abs(rh[peaks,1])))[0][0]
    for i in range(mx,len(peaks)):
        if abs(rh[peaks[i],1])<abs(rh[peaks[i+1],1]):
            mn=i
            break
    for i in range(mn,len(peaks)):
        if abs(rh[peaks[i],1])>abs(rh[peaks[i+1],1]):
            mx2=i
            break
    dt=rh[peaks[mn],0]-rh[peaks[mx],0]
    rhM=abs(rh[peaks[mx2],1])
    return dt,rhM


#BAM data analysis
dtM2=np.zeros(31)
rhM2=np.zeros(31)
names2=list()
code2=np.zeros(31)
l=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'data/BAM:0%s%s%s.h5' %(m,k,j)
            try:
                file=h5py.File(name,'r')
                dat = list(file["/rh_22"])
                rh = np.array(file["/rh_22/%s" %dat[-1]])
                dtM2[l],rhM2[l]=analyze(rh)
                names2.append(name)
                code2[l]=m*100+k*10+j
                l=l+1

            except OSError:
                pass

#THC data analysis
dtM1=np.zeros(23)
rhM1=np.zeros(23)
names1=list()
code1=np.zeros(23)
l=0
for m in range(0,4):
    for k in range(0,10):
        for j in range(0,10):
            name = 'data/THC:0%s%s%s.h5' %(m,k,j)
            try:
                file=h5py.File(name,'r')
                dat = list(file["/rh_22"])
                rh = np.array(file["/rh_22/Rh_l2_m2_r00400.txt" ])
                dtM1[l],rhM1[l]=analyze(rh)
                names1.append(name)
                code1[l]=m*100+k*10+j
                l=l+1

            except OSError:
                pass
dtM=np.concatenate((dtM1,dtM2))
rhM=np.concatenate((rhM1,rhM2))
names=np.concatenate((names1,names2))
code=np.concatenate((code1,code2))


#BAM metadata

SLy=1
H4=2
MS1=3
MPA1=4
ALF2=5
MS1b=6
ms1b=6
ENG=7

q2 = np.zeros(31)
mas2=np.zeros(31)
eos=np.zeros(31)


i=0
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'metadata/BAM:0%s%s%s.txt' %(m,k,j)
            try:
                f=open(name)
                lines=f.readlines()
                exec(lines[10])
                q2[i] = id_mass_ratio
                exec(lines[8])
                mas2[i]=id_mass
                if i==0:
                    eos[i]=8 #by hand because 2H can not be given a value
                if i>0:
                    exec(lines[15])
                    eos[i]=id_eos

                i=i+1
            except OSError:
                pass
mas12=np.zeros([31,2])
for i in range(31):
    mas12[i,0]=(q2[i]*mas2[i])/(1+q2[i])
    mas12[i,1]=(mas2[i])/(1+q2[i])

#Tidal Deformability
m_r1=np.load('tid_def/SLy.npy')
m_r2=np.load('tid_def/H4.npy')
m_r3=np.load('tid_def/MS1.npy')
m_r4=np.load('tid_def/MPA1.npy')
m_r5=np.load('tid_def/ALF2.npy')
m_r6=np.load('tid_def/MS1b.npy')
m_r7=np.load('tid_def/ENG.npy')
m_r8=np.load('tid_def/2H.npy')


k_l1=np.load('tid_def/k_l_SLy.npy')
k_l2=np.load('tid_def/k_l_H4.npy')
k_l3=np.load('tid_def/k_l_MS1.npy')
k_l4=np.load('tid_def/k_l_MPA1.npy')
k_l5=np.load('tid_def/k_l_ALF2.npy')
k_l6=np.load('tid_def/k_l_MS1b.npy')
k_l7=np.load('tid_def/k_l_ENG.npy')
k_l8=np.load('tid_def/k_l_2H.npy')



mx=np.amax(m_r1[0])
idx=np.where(m_r1[0]==mx)
idx=idx[0][0]
cs1=spline(m_r1[0][1:idx],k_l1[0][1:idx])
cs11=spline(m_r1[0][1:idx],m_r1[1][1:idx])

mx=np.amax(m_r2[0])
idx=np.where(m_r2[0]==mx)
idx=idx[0][0]
cs2=spline(m_r2[0][1:idx],k_l2[0][1:idx])
cs21=spline(m_r2[0][1:idx],m_r2[1][1:idx])

mx=np.amax(m_r3[0])
idx=np.where(m_r3[0]==mx)
idx=idx[0][0]
cs3=spline(m_r3[0][1:idx],k_l3[0][1:idx])
cs31=spline(m_r3[0][1:idx],m_r3[1][1:idx])

mx=np.amax(m_r4[0])
idx=np.where(m_r4[0]==mx)
idx=idx[0][0]
cs4=spline(m_r4[0][1:idx],k_l4[0][1:idx])
cs41=spline(m_r4[0][1:idx],m_r4[1][1:idx])

mx=np.amax(m_r5[0])
idx=np.where(m_r5[0]==mx)
idx=idx[0][0]
cs5=spline(m_r5[0][1:idx],k_l5[0][1:idx])
cs51=spline(m_r5[0][1:idx],m_r5[1][1:idx])

mx=np.amax(m_r6[0])
idx=np.where(m_r6[0]==mx)
idx=idx[0][0]
cs6=spline(m_r6[0][1:idx],k_l6[0][1:idx])
cs61=spline(m_r6[0][1:idx],m_r6[1][1:idx])

mx=np.amax(m_r7[0])
idx=np.where(m_r7[0]==mx)
idx=idx[0][0]
cs7=spline(m_r7[0][1:idx],k_l7[0][1:idx])
cs71=spline(m_r7[0][1:idx],m_r7[1][1:idx])

mx=np.amax(m_r8[0])
idx=np.where(m_r8[0]==mx)
idx=idx[0][0]
cs8=spline(m_r8[0][1:idx],k_l8[0][1:idx])
cs81=spline(m_r8[0][1:idx],m_r8[1][1:idx])



k212=np.zeros([31,2])
r12=np.zeros([31,2])
for i in range(31):
    if eos[i]==1:
        k212[i,0]=cs1(mas12[i,0])
        k212[i,1]=cs1(mas12[i,1])
        r12[i,0]=cs11(mas12[i,0])
        r12[i,1]=cs11(mas12[i,1])
    elif eos[i]==2:
        k212[i,0]=cs2(mas12[i,0])
        k212[i,1]=cs2(mas12[i,1])
        r12[i,0]=cs21(mas12[i,0])
        r12[i,1]=cs21(mas12[i,1])
    elif eos[i]==3:
        k212[i,0]=cs3(mas12[i,0])
        k212[i,1]=cs3(mas12[i,1])
        r12[i,0]=cs31(mas12[i,0])
        r12[i,1]=cs31(mas12[i,1])
    elif eos[i]==4:
        k212[i,0]=cs4(mas12[i,0])
        k212[i,1]=cs4(mas12[i,1])
        r12[i,0]=cs41(mas12[i,0])
        r12[i,1]=cs41(mas12[i,1])
    elif eos[i]==5:
        k212[i,0]=cs5(mas12[i,0])
        k212[i,1]=cs5(mas12[i,1])
        r12[i,0]=cs51(mas12[i,0])
        r12[i,1]=cs51(mas12[i,1])
    elif eos[i]==6:
        k212[i,0]=cs6(mas12[i,0])
        k212[i,1]=cs6(mas12[i,1])
        r12[i,0]=cs61(mas12[i,0])
        r12[i,1]=cs61(mas12[i,1])
    elif eos[i]==7:
        k212[i,0]=cs7(mas12[i,0])
        k212[i,1]=cs7(mas12[i,1])
        r12[i,0]=cs71(mas12[i,0])
        r12[i,1]=cs71(mas12[i,1])
    elif eos[i]==8:
        k212[i,0]=cs8(mas12[i,0])
        k212[i,1]=cs8(mas12[i,1])
        r12[i,0]=cs81(mas12[i,0])
        r12[i,1]=cs81(mas12[i,1])


l=np.zeros([31,2])
lamda2=np.zeros(31)
for i in range(31):
    l[i,0]=(2.0/3.0)*k212[i,0]*(1.0/pow(mas12[i,0]/r12[i,0],5))
    l[i,1]=(2.0/3.0)*k212[i,1]*(1.0/pow(mas12[i,1]/r12[i,1],5))
    lamda2[i]=(16/13)*( ( (mas12[i,0]+12*mas12[i,1])*pow(mas12[i,0],4)*l[i,0]+(mas12[i,1]+12*mas12[i,0])*pow(mas12[i,1],4)*l[i,1] )/pow( mas12[i,0]+mas12[i,1],5 ) )


#THC metadata
q1 = np.zeros(23)
mas1=np.zeros(23)
lamda1=np.zeros(23)
i=0
for m in range(0,4):
    for k in range(0,10):
        for j in range(0,10):
            name = 'metadata/THC:0%s%s%s.txt' %(m,k,j)
            try:
                f=open(name)
                lines=f.readlines()
                exec(lines[10])
                q1[i] = id_mass_ratio
                exec(lines[8])
                mas1[i]=id_mass
                exec(lines[17])
                lamda1[i]=id_Lambda
                i=i+1
            except OSError:
                pass


lamda=np.concatenate((lamda1,lamda2))
q=np.concatenate((q1,q2))
mas=np.concatenate((mas1,mas2))


if os.path.exists('results'):
    pass
else:
    os.mkdir('results')

np.save('results/lamda.npy',lamda)
np.save('results/q.npy',q)
np.save('results/rhM.npy',rhM)
np.save('results/dtM.npy',dtM)

fig1=plt.figure()
plt.scatter(lamda,dtM)
plt.xlabel('lamda')
plt.ylabel('dtM')
plt.show()

fig2=plt.figure()
plt.scatter(q,rhM)
plt.xlabel('q')
plt.ylabel('rhM')
plt.show()
