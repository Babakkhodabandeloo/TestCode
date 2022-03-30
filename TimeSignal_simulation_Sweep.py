import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.fftpack
import json
#PI=math.pi
#print(2*PI)


#%%
SaveJSON_Dir='/home/bkh/OnlineCourses/BroadBandBackScatt_TimeSignal/SaveJSON_Dir/'

import os
os.chdir('/home/bkh/OnlineCourses/PYfunctions/')
from FUNC_generate_SweepSig import func_gen_sweepsig

f0=70000
Deltaf=20000
T=2*1.024
mexpo=20
[t_vec,p_vec]=func_gen_sweepsig(f0,Deltaf,T,mexpo)

if np.mod(len(t_vec),2)==1:
   t_vec=t_vec[0:len(t_vec)-1]  
p_vec=p_vec[0:len(t_vec)]

FIG=plt.figure(figsize=(6,4))  
X_shifts=[0.06]
Y_shifts=[0.03]
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
ax = FIG.add_subplot(111)
box = ax.get_position()
#        ax.set_position([box.x0+0.0, box.y0+(0.01), box.width * 1 , box.height * 1.0])
ax.set_position([box.x0+X_shifts[0], box.y0+Y_shifts[0], box.width * 1.0 , box.height * 1.0])
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
plt.plot(1E3*t_vec,p_vec,color=[0,0,1])
plt.xticks(fontsize=14, rotation=0)
plt.yticks(fontsize=14, rotation=0)
plt.xlabel('time (ms)',fontsize=14)
plt.ylabel('Pressure (kPa) ',fontsize=14) 
#plt.show()


#%% Fourier Transform
dt=np.mean(t_vec[1:]-t_vec[0:-1])
Fs=1/dt

DF=1/(dt*len(t_vec))
F_vec = np.fft.fftfreq(len(t_vec),1)/(dt)
#plt.plot(F_vec)
#F_vec=np.arange(0,Fs+DF,DF)
F_vec=F_vec[0:len(t_vec)]
#plt.plot(F_vec)
DFT_p=scipy.fftpack.fft(p_vec)
FIG=plt.figure(figsize=(6,4))  
plt.plot(scipy.fft.fftshift(F_vec)/1000,scipy.fft.fftshift(abs(DFT_p)),'k')


FIG=plt.figure(figsize=(6,4))  
plt.plot(scipy.fft.fftshift(F_vec)/1000,scipy.fft.fftshift(abs(DFT_p))/np.max(abs(DFT_p)),'k')
plt.xlabel('Frequency (kHz)')
plt.title('Normalized Fourier Transform of Transmitt Signal')

IDFT=scipy.fftpack.ifft(DFT_p)
TransmitSig=np.real(IDFT)
#%% Backscattering
import os
os.chdir('/home/bkh/OnlineCourses/PYfunctions/')
from FUNC_fromfunc_ElasticSphere_WC import func_formfunc_ElasticSphere_WC
#func_BackScatt_Sphere_Anderson_f_bs(ro_water,ro_bubble,c_water,c_bubble,bubble_radius_m,freq_kHz)

LocVec=np.where( (abs(DFT_p[0:int(0.5*len(DFT_p))])/np.max(abs(DFT_p)) ) >0.05)[0]
Ind1=LocVec[0]
Ind2=LocVec[-1]
#Ind1=0
#Ind2=int(0.5*len(DFT_p))-5

DFT_p_bs=0*DFT_p # this will be modified in the loop coming next
DFT_reconst=0*DFT_p

TS_vec=[]
Freqvec=[]

#ro_water=1000
#ro_bubble=80#1000*1.05
#c_water=1500
#c_bubble=325 #1500*1.05
#bubble_radius_m=0.00065

M_order=int(0.5*38.1*1E-3*2*np.pi*F_vec[Ind2]/1500)+20
R=0

for ii in range(Ind1,Ind2+1):
#    print(F_vec[ii])
    Freqvec=np.append(Freqvec,F_vec[ii])
    F_bs=func_formfunc_ElasticSphere_WC(F_vec[ii],M_order)
#    F_bs=0.1-0.15j
#    F_bs=0.2-0.1j
#    DFT_p_bs[ii]=DFT_p[ii]*F_bs
#    DFT_p_bs[int(len(DFT_p_bs))-ii-1]=DFT_p[int(len(DFT_p_bs))-ii-1]*np.conj(F_bs)
#    TS_vec=np.append(TS_vec,20*np.log10(abs(F_bs)))
#    
#    DFT_reconst[ii]=DFT_p_bs[ii]*(1/F_bs)
#    DFT_reconst[int(len(DFT_p_bs))-ii-1]=DFT_p_bs[int(len(DFT_p_bs))-ii-1]*np.conj(1/F_bs)
    w=2*np.pi*F_vec[ii]
    DFT_p_bs[ii]=DFT_p[ii]*np.exp(-1j*w*R/1500)*F_bs
    DFT_p_bs[int(len(DFT_p_bs))-ii]=DFT_p[int(len(DFT_p_bs))-ii]*np.exp(-1j*(-w)*R/1500)*np.conj(F_bs)
    TS_vec=np.append(TS_vec,20*np.log10(abs(F_bs)))
    
    DFT_reconst[ii]=DFT_p_bs[ii]*np.exp(1j*w*R/1500)*(1/F_bs)
    DFT_reconst[int(len(DFT_p_bs))-ii]=DFT_p_bs[int(len(DFT_p_bs))-ii]*np.exp(1j*(-w)*R/1500)*np.conj(1/F_bs)
    
#    DFT_p_bs[ii]=DFT_p[ii]*F_bs
#    DFT_p_bs[int(len(DFT_p_bs))-ii]=DFT_p[int(len(DFT_p_bs))-ii]*np.conj(F_bs)
    
##plt.plot(abs(DFT_p_bs))  
    
FIG=plt.figure(figsize=(8,7)) 
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
ax = FIG.add_subplot(411)
box = ax.get_position()
#        ax.set_position([box.x0+0.0, box.y0+(0.01), box.width * 1 , box.height * 1.0])
ax.set_position([box.x0+X_shifts[0], box.y0+Y_shifts[0], box.width * 1.0 , box.height * 1.0])
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
plt.plot(1E3*t_vec,np.real(IDFT),color='r',dashes=[2,0],linewidth=1.5)
    
IDFT_p_bs=scipy.fftpack.ifft(DFT_p_bs)
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
ax = FIG.add_subplot(412)
box = ax.get_position()
#        ax.set_position([box.x0+0.0, box.y0+(0.01), box.width * 1 , box.height * 1.0])
ax.set_position([box.x0+X_shifts[0], box.y0+Y_shifts[0], box.width * 1.0 , box.height * 1.0])
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
BackScattSignal=np.real(IDFT_p_bs)
plt.plot(1E3*t_vec,BackScattSignal,color='k',dashes=[2,0],linewidth=1.5)

#plt.plot(abs(DFT_p_bs))    
IDFT_reconst=scipy.fftpack.ifft(DFT_reconst)
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
ax = FIG.add_subplot(413)
box = ax.get_position()
#        ax.set_position([box.x0+0.0, box.y0+(0.01), box.width * 1 , box.height * 1.0])
ax.set_position([box.x0+X_shifts[0], box.y0+Y_shifts[0], box.width * 1.0 , box.height * 1.0])
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
plt.plot(1E3*t_vec,np.real(IDFT_reconst),color=[0.5,0.5,0.5],dashes=[2,0],linewidth=1.5)

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  
ax = FIG.add_subplot(414)
box = ax.get_position()
#        ax.set_position([box.x0+0.0, box.y0+(0.01), box.width * 1 , box.height * 1.0])
ax.set_position([box.x0+X_shifts[0], box.y0+Y_shifts[0], box.width * 1.0 , box.height * 1.0])
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 
plt.plot(Freqvec,TS_vec,color=[1,0.2,0.2])
plt.show()

#FIG=plt.figure(figsize=(6,4))  
#plt.plot(scipy.fft.fftshift(F_vec)/1000,scipy.fft.fftshift(abs(DFT_p))/np.max(abs(DFT_p)),'k')
#plt.plot(scipy.fft.fftshift(F_vec)/1000,scipy.fft.fftshift(abs(DFT_p_bs))/np.max(abs(DFT_p_bs)),color='r',dashes=[4,2],linewidth=1.5)

#%% Preparing Signals
FIG=plt.figure(figsize=(8,7)) 
M_order=int(0.5*38.1*1E-3*2*np.pi*F_vec[Ind2]/1500)+20
R=0

Freqvec=[]
TS_vec=[]
DFT_p_bs=0*DFT_p
for ii in range(Ind1,Ind2+1):
#    print(F_vec[ii])
    Freqvec=np.append(Freqvec,F_vec[ii])
    F_bs=func_formfunc_ElasticSphere_WC(F_vec[ii],M_order)
#    F_bs=0.1-0.15j
#    F_bs=0.2-0.1j
#    DFT_p_bs[ii]=DFT_p[ii]*F_bs
#    DFT_p_bs[int(len(DFT_p_bs))-ii-1]=DFT_p[int(len(DFT_p_bs))-ii-1]*np.conj(F_bs)
#    TS_vec=np.append(TS_vec,20*np.log10(abs(F_bs)))
#    
#    DFT_reconst[ii]=DFT_p_bs[ii]*(1/F_bs)
#    DFT_reconst[int(len(DFT_p_bs))-ii-1]=DFT_p_bs[int(len(DFT_p_bs))-ii-1]*np.conj(1/F_bs)
    w=2*np.pi*F_vec[ii]
    DFT_p_bs[ii]=DFT_p[ii]*np.exp(-1j*w*R/1500)*F_bs
    DFT_p_bs[int(len(DFT_p_bs))-ii]=DFT_p[int(len(DFT_p_bs))-ii]*np.exp(-1j*(-w)*R/1500)*np.conj(F_bs)
    TS_vec=np.append(TS_vec,20*np.log10(abs(F_bs)))
    
IDFT_p_bs=scipy.fftpack.ifft(DFT_p_bs)
BackScattSignal1=np.real(IDFT_p_bs)
plt.plot(1E3*t_vec,BackScattSignal1,color='k',dashes=[2,0],linewidth=1.5)

M_order=int(0.5*38.1*1E-3*2*np.pi*F_vec[Ind2]/1500)+20
R=0.25

Freqvec=[]
TS_vec=[]
DFT_p_bs=0*DFT_p
for ii in range(Ind1,Ind2+1):
#    print(F_vec[ii])
    Freqvec=np.append(Freqvec,F_vec[ii])
    F_bs=func_formfunc_ElasticSphere_WC(F_vec[ii],M_order)
#    F_bs=0.1-0.15j
#    F_bs=0.2-0.1j
#    DFT_p_bs[ii]=DFT_p[ii]*F_bs
#    DFT_p_bs[int(len(DFT_p_bs))-ii-1]=DFT_p[int(len(DFT_p_bs))-ii-1]*np.conj(F_bs)
#    TS_vec=np.append(TS_vec,20*np.log10(abs(F_bs)))
#    
#    DFT_reconst[ii]=DFT_p_bs[ii]*(1/F_bs)
#    DFT_reconst[int(len(DFT_p_bs))-ii-1]=DFT_p_bs[int(len(DFT_p_bs))-ii-1]*np.conj(1/F_bs)
    w=2*np.pi*F_vec[ii]
    DFT_p_bs[ii]=DFT_p[ii]*np.exp(-1j*w*R/1500)*F_bs
    DFT_p_bs[int(len(DFT_p_bs))-ii]=DFT_p[int(len(DFT_p_bs))-ii]*np.exp(-1j*(-w)*R/1500)*np.conj(F_bs)
    TS_vec=np.append(TS_vec,20*np.log10(abs(F_bs)))

#ii=501
#print(DFT_p[ii])
#print(DFT_p[int(len(DFT_p_bs))-ii])

    
IDFT_p_bs=scipy.fftpack.ifft(DFT_p_bs)
BackScattSignal2=np.real(IDFT_p_bs)
plt.plot(1E3*t_vec,BackScattSignal2,color='b',dashes=[2,0],linewidth=1.5)

BackScattSignal=BackScattSignal1+0*BackScattSignal2
plt.plot(1E3*t_vec,BackScattSignal,color='k',dashes=[2,0],linewidth=1.5)
#%% Write JSON file
DB={}
cnt=0
DB[cnt]={}
#DB[cnt]['TS_Num']=TS_Num_Vec[cnt]
DB[cnt]['Freqvec']=Freqvec.tolist()
DB[cnt]['TS_vec']=TS_vec.tolist()
DB[cnt]['TimeSignal']=t_vec.tolist()
DB[cnt]['BackScattSignal']=BackScattSignal.tolist()
DB[cnt]['TransmittSignal']=TransmitSig.tolist()


JsonOutfile='TimeSignalSimulation_Distance0'+str(R)+'m.json'
#JsonOutfile='TimeSignalSimulation_SingleTarget.json'
SaveJsonDirFile=SaveJSON_Dir+JsonOutfile
with open(SaveJsonDirFile, 'w') as Outfile:
    json.dump(DB,Outfile,sort_keys=True, indent=2) 
#%% Read JSON
#list(DB.keys())
#list(DB[cnt].keys())    
#DB[cnt]['Eta2'] 
    
    
#%% FFT of backscattered and transmit signal
FIG=plt.figure(figsize=(8,7)) 
DFT_BackScattSignal=scipy.fftpack.fft(BackScattSignal)
ax = FIG.add_subplot(311)
plt.plot(Freqvec[:-1]/1000,10*np.log10(np.abs(DFT_BackScattSignal[Ind1:Ind2])**2))

DFT_TransmitSig=scipy.fftpack.fft(TransmitSig)
ax = FIG.add_subplot(312)
plt.plot(Freqvec[:-1]/1000,10*np.log10(np.abs(DFT_TransmitSig[Ind1:Ind2])**2))

Ratio=np.abs(DFT_BackScattSignal[Ind1:Ind2])/np.abs(DFT_TransmitSig[Ind1:Ind2])

ax = FIG.add_subplot(313)
plt.plot(Freqvec[:-1]/1000,10*np.log10(np.abs(Ratio)**2),color=[0,0,1],linewidth=1.5)
plt.plot(Freqvec/1000,TS_vec,color=[1,0.2,0.2],dashes=[3,3],linewidth=2)