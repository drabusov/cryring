from matplotlib import *
import os
import matplotlib.pyplot as plt
import math
import numpy as np
import random
from numpy.fft import rfft

from scipy.interpolate import interp1d

#---------------------------------------------PLOTTING-SETTINGS--------------------------------------
rc('text', usetex=False)
rc('font', serif ='Times')
rc('font', size=16)
rc('axes', linewidth=0.5)
rc('lines', linewidth=1.15	)
#rc('figure', figsize=(8.3,5.2))

#rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
##Options
#params = {'text.usetex' : True,
#          'font.size' : 18,
#          'font.family' : 'lmodern',
#          'text.latex.unicode': True,
#          'figure.figsize' : (6,5)
#          }
#rcParams.update(params) 
#----------------------------------------------------------------------------------------------------

#--------------------------------------------SAVING-PICTURES---------------------------------------
def save(name=''):
    pwd = os.getcwd()
    os.chdir('./Pictures')                    # Folder ./Pictures must exist
    plt.savefig('%s.pdf' % name, fmt ='pdf')
    os.chdir(pwd)
#----------------------------------------------------------------------------------------------------




SP2 = 48


def read_twiss(filename, header):
    
    twiss = open(filename)  # This file updates every time of madx running

    s, sbpm, beta, betabpm, phasebpm = list(), list(), list(), list(), list()

    for i,line in enumerate(twiss):
        line=line.strip()
        line=line.split()
        if line[1] == 'Q1':
            tune = float(line[3])
        if i > header:
            s.append(float(line[2]))
            beta.append(float(line[6])) 
            if 'DX1' in line[0]:
                sbpm.append(float(line[2]))
                betabpm.append(float(line[6]))
                phasebpm.append(2*np.pi*float(line[10]))
# what?           dr = float(line[2]) - svar[i-SP2-3]

    print(i)
    twiss.close()
    return tune, s, sbpm, beta, betabpm, phasebpm

Q1, svar, sbpm, betamod, betabpm, phasebpm = read_twiss('twiss.txt',SP2)



ss = [svar[0]]
bb = [betamod[0]]
i = 1
while i < len(svar):
    if svar[i] !=svar[i-1]:
        ss.append(svar[i])
        bb.append(betamod[i])
    i+=1


s_inter = np.linspace(svar[0], svar[-1], len(svar)*2)
inter = interp1d(ss, bb, 'cubic')
beta_inter = inter(s_inter)

#plt.figure()
#plt.scatter(svar,betamod, s = 10)
#plt.plot(s_inter,beta_inter)




beta_spec = np.abs(rfft(beta_inter-np.mean(beta_inter)))
qq = np.linspace(0,2/(s_inter[1]-s_inter[0]), len(s_inter))
frq = qq[6]
beta_sin = -2*beta_spec[6]/len(s_inter)*np.cos(2*np.pi*frq*s_inter) +np.mean(beta_inter)

#fig, axs = plt.subplots(1, 2)
#fig.subplots_adjust(hspace=0.4,wspace = 0.3)
#axs[0].plot(s_inter,beta_inter)
#axs[0].plot(s_inter,beta_sin)
#axs[1].plot(beta_spec)
#axs[1].grid(True)



print('horizontal tune (from TWISS) is {}'.format(Q1))



dr = list()

phase = list()
phase1 = list()

phaserror = list()
betaerror = list()

phaserror2 = list()
betaerror2 = list()

phaserror12 = list()
betaerror12 = list()

phaserror22 = list()
betaerror22 = list()


#-----------------------------------TRACK-READING-----------------------------------------------
f = open('record.txtone')

M = len(betabpm)         # Number of BPMs
SP = 54          # lines without info of track in trackone file 

i=0
k=0

x,y,tau,pt = list(),list(), list(), list()
for i,line in enumerate(f):
    if i > SP and i%2!=1:
        k = k+1
        if k%(M+1) != 0:
            line=line.strip()
            line=line.split()
            x.append(float(line[2]))
            y.append(float(line[4]))
            tau.append(float(line[6]))
            pt.append(float(line[7]))


N = k//(M+1) # Number of turns
A = max(x) # Amplitude. If the orbit is zero
print('Number of turns = {}, number of BPMs = {}'.format(N,M))

phasemod = list()
for i in range(1,M):
    phasemod.append(phasebpm[i]-phasebpm[i-1])
#    phasemod2.append(phasebpm2[i]-phasebpm2[i-1])

x = np.reshape(x, (N,M))
m = np.transpose(x)

y = np.reshape(y, (N,M))
my = np.transpose(y)          # For FFT (cheking vertical tune)

tau = np.reshape(tau, (N,M))
mt = np.transpose(tau)          # For FFT (cheking vertical tune)

pt = np.reshape(pt, (N,M))
mp = np.transpose(pt)          # For FFT (cheking vertical tune)


plt.figure()
#plt.plot(mt[0]/max(mt[0]))
plt.plot(mt[0])
plt.plot(mp[0]/min(mp[0]))
plt.grid(True)

#--------------------------------------SINGLE-RUN-------------------------------------------

#noise = list()
#for i in range(0,N*M):
#    noise.append(random.uniform(-1,1))

#noise = np.reshape(noise,(M,N))

#m2 = range(M)

#m1 = m + 0.1*noise               # noise = 100*10^-6 meters


m1 = [x-np.mean(x) for x in m]

#for i in range(M):
#    m1[i] = m1[i]-np.mean(m1[i]) 
#    m2[i] = m1[i][0:200]         # Cut the length of beam history (200 turns)


U, s, V = np.linalg.svd(np.transpose(m1), full_matrices=True)

u1 = np.transpose(U)
mode = s[0]*u1[0]+s[1]*u1[1]

spec_mode = np.abs(rfft(mode))
q = np.linspace(0,0.5,len(spec_mode))


plt.figure()
plt.plot(mode)
plt.grid(True)
plt.xlabel("n, turn index")
plt.title("Betatron eigenmode. MADX")
#plt.savefig("mode_ripple.pdf")


plt.figure()
plt.plot(q,spec_mode)
plt.xlim([0.41,0.43])
plt.grid(True)
plt.xlabel("tune")
plt.title("Spectrum of betatron eigenmode. MADX")
#plt.savefig("spec_mode_ripple.pdf")


fig, axs = plt.subplots(4, 2)
fig.subplots_adjust(hspace=0.4,wspace = 0.3)
[axs[i//2,i%2].plot(x) for i,x in enumerate(m)]

beta = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)
scale = betabpm[0]/beta[0]   # put the correct scaling factor by hand
beta = beta*scale

for i in range(0,M):
    phase.append(math.atan(s[1]*V[1][i]/(s[0]*V[0][i])))

for k in range(1,M):
    if abs(phase[k]-phase[k-1]) > np.pi/2.:
        phase1.append(np.pi - abs(phase[k]-phase[k-1]))
    else:
        phase1.append(abs(phase[k]-phase[k-1])) 


#-------------------------------FFT--------------------------------------------------------------

y = rfft(m[1])
yz = rfft(my[1])
p = np.argmax(y)
py = np.argmax(yz)
Q = np.linspace(0,0.5,len(y))
print('FFT. horizontal tune is {} OR {}'.format(Q[p],1-Q[p]))
print('FFT. vertical tune is {} OR {}'.format(Q[py],1-Q[py]))        # tune and mirrored tune

#------------------------------------------------------------------------------------------------

#-----------------------------------NOISE-OF-BPM-------------------------------------------------
'''
noise = list()

#for i in range(0,N*M):
#    noise.append(random.uniform(-1,1))         # constructing a noise (sequence of noncorrelated numbers) for 2 (!!) cases

#for i in noise:
#    noi.write(repr(i)+'\n')                    # if it works next loop must be commented

for line in noi:
    line = line.strip()
    noise.append(float(line))                  # if it works 2 (!!) previous loops must be commented

noise = np.reshape(noise,(M,N))

per = np.array(range(0, 300, 5))

t=0
for q in per:
    phase = list()
    phase1 = list()
    m2 = range(M)
    m1 = m + q*0.001*noise
    for i in range(0,M):
        m1[i] = m1[i]-np.mean(m1[i])
        m2[i] = m1[i][0:200]
    U, s, V = np.linalg.svd(np.transpose(m2), full_matrices=True)
    beta = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)
    beta = beta*betabpm[0]/beta[0]
    betaerror.append(100*np.mean(abs(np.array(beta)- np.array(betabpm))/np.array(betabpm)))
    for i in range(0,M):
        phase.append(math.atan(s[1]*V[1][i]/(s[0]*V[0][i])))
    for k in range(1,M):
#        phase1.append(np.sign(phase[k]-phase[k-1]))
        if abs(phase[k]-phase[k-1]) > np.pi/2.:
            phase1.append(np.pi - abs(phase[k]-phase[k-1]))
        else:
            phase1.append(abs(phase[k]-phase[k-1]))
    phaserror.append(100*np.mean(abs(np.array(phase1)-np.array(phasemod))/np.array(phasemod)))
#    fn2.write(repr(q)+' '+repr(betaerror[t])+' '+repr(phaserror[t])+' '+repr(Q1)+'\n')
#    t = t+1
'''
#----------------------------------------------------------------------------------------------

#------------------------------NUMBER-OF-TURNS-------------------------------------------------
'''
noise = list()
m2 = list()

per = np.array(range(10,250))

#for i in range(0,N*M):
#    noise.append(random.uniform(-1,1))         # constructing a noise (sequence of noncorrelated numbers) for 2 (!!) cases

#for i in noise:
#    noi.write(repr(i)+'\n')                    # if it works next loop must be commented

for line in noi:
    line = line.strip()
    noise.append(float(line))                  # if it works previous loop must be commented

noise = np.reshape(noise,(M,N))

t=0
for P in per:
    phase = list()
    phase1 = list()
    m2 = range(M)
    m1 = m + 0.1*noise
    for i in range(0,M):
        m1[i] = m1[i]-np.mean(m1[i])
        m2[i] = m1[i][0:P]
    U, s, V = np.linalg.svd(np.transpose(m2), full_matrices=True)
    beta = pow(s[0],2)*(V[0]**2)+pow(s[1],2)*(V[1]**2)
    beta = beta*betabpm[11]/beta[11]
    betaerror.append(100*np.mean(abs(np.array(beta)- np.array(betabpm))/np.array(betabpm)))
    for i in range(0,M):
        phase.append(math.atan(s[1]*V[1][i]/(s[0]*V[0][i])))
    for k in range(1,M):
#        phase1.append(np.sign(phase[k]-phase[k-1]))
        if abs(phase[k]-phase[k-1]) > np.pi/2.:
            phase1.append(np.pi - abs(phase[k]-phase[k-1]))
        else:
            phase1.append(abs(phase[k]-phase[k-1]))
    phaserror.append(100*np.mean(abs(np.array(phase1)-np.array(phasemod))/np.array(phasemod))) 
    ft2.write(repr(P)+' '+repr(betaerror[t])+' '+repr(phaserror[t])+'\n')
    t = t+1
'''
#---------------------------------------------------------------------------------------------------

#----------------------------------PLOTTING-PART----------------------------------------------------

#-------------------------------PLOT---------BETA---------------------------------------------------
'''
fig = plt.figure()

gr6 = plt.scatter(sbpm,betabpm2,color='b',label = r'Ideal $\beta$-function (MADX)')#r'MADX $\beta$-function')
gr6 = plt.scatter(sbpm,betabpm,color='g', marker = '^',s=50,label = r'Distorted $\beta$-function (MADX)')
gr5 = plt.scatter(sbpm,beta,color='r',marker = 'x', s=75, label = r'MIA $\beta$-function')
text1 = plt.text(100,21.5,r'$\nu_x = 18.84$'+'\n'+r'$\nu_y = 18.73$'+'\n'+r'$\Delta x/A = 0.05$')
text2 = plt.text(30.6,38,'correct scaling factor = %.3f'%scale+'\n'+'gradient error = 5%, 8 quads')

plt.grid(True)
#axes = plt.gca()
#axes.set_xlim([-15,1100])
#axes.set_ylim([8,45])
plt.legend(loc='upper right',frameon=True)

plt.title(r'$\beta$-function at BPMs of SIS100')
plt.xlabel('s, m')
plt.ylabel(r'$\beta$, m')

#save(name='beta')
'''
#------------------------------PLOT-----SECTOR------------------------------------------------------------------

fig = plt.figure()

plt.plot(svar,betamod,color = 'black',label = r'Ideal $\beta$ (MADX)')
#plt.plot(svar,betamod, color = 'g', ls ='--', label =r'Distorted $\beta$ (MADX)') #label = r'MADX $\beta$-function')#
#gr6 = plt.scatter(sbpm,betabpm2,color='b',label = r'Ideal on BPMs')
gr6 = plt.scatter(sbpm,betabpm,color='g',s=50, marker = '^',label = r'MADX $\beta$ at BPMs' )#,label = r'Distorted on BPMs')
gr5 = plt.scatter(sbpm,beta,color='r',marker = 'x', s=95, label = r'MIA $\beta$-function')

#plt.text(180.6*2+5,37,'8 quadrupoles: 5% error of gradient')
#text2 = plt.text(180.6*2+5,35,'incorrect scaling factor = %.3f'%scale)

#text1 = plt.text(180.6*2+15,23,r'$\nu_x = 18.84$'+'\n'+r'$\nu_y = 18.73$'+'\n'+r'$\Delta x/A = 0.05$')

plt.grid(True)
axes = plt.gca()
axes.set_xlim([svar[0],svar[-1]])
axes.set_ylim([min(betamod) - max(betamod)/15,max(betamod)+max(betamod)/2])
plt.legend(loc='upper right',frameon=True)

plt.title(r'$\beta$-function (one sector)')
plt.xlabel('s, m')
plt.ylabel(r'$\beta$, m')

#save(name='beta_sector')

#------------------------------PLOT---------PHASE---------------------------------------------------------------
'''
fig = plt.figure()

gr9 = plt.scatter(range(2,M+1),phasemod, color = 'b', s=10, label = r'MADX phase advance')
gr8 = plt.scatter(range(2,M+1),phase1, color = 'r',marker = '+', s=80, label = r'MIA phase advance')
text1 = plt.text(7,1.6,r'$\nu_x = 18.84$'+'\n'+r'$\nu_y = 18.73$'+'\n'+r'$\Delta x/A = 0.05$')

plt.grid(True)
axes = plt.gca()
axes.set_xlim([-5,90])
axes.set_ylim([1.15,1.73])
plt.legend(loc='best',frameon=True)

plt.title(r'Phase advance',fontsize=18)
plt.xlabel('BPM, i',fontsize=18)
plt.ylabel(r'$\Psi_i - \Psi_{i-1}$, rad',fontsize=18)

#save(name='phase_present_3')
'''
#------------------------PLOT--------------TRACKS-------------------------------------------------------------------
'''
fig = plt.figure()
time = np.array(range(0,50))
gr1 = plt.plot(m[0][50:100],color = 'r', label = 'BPM-1')
gr2 = plt.plot(m[4][50:100],color = 'b', label = 'BPM-5')
#gr3 = plt.plot(m[8][50:100], label = 'BPM-2')

ampl = 2.
coef = np.sqrt(betabpm[0]/betabpm[8])
#coefm = np.sqrt(beta[0]/beta[8])

#plt.plot(range(0,50,49),[ ampl*coefm, ampl*coefm],color='r',ls ='--')
#plt.plot(range(0,50,49),[-ampl*coefm,-ampl*coefm],color='r',ls ='--')

plt.plot(range(0,50,49),[ ampl*coef, ampl*coef],color='black',ls ='--')
plt.plot(range(0,50,49),[-ampl*coef,-ampl*coef],color='black',ls ='--')

#plt.annotate(r'$A =$'+repr(ampl)+r'$\cdot \sqrt{\frac{\beta 1}{\beta 9}}$',fontsize = 20, xy=(30, ampl*coef), xytext=(2, 2.3), arrowprops=dict(facecolor='black', shrink=0.05),)

plt.plot(range(0,50,49),[ ampl, ampl],color='b',ls ='--')
plt.plot(range(0,50,49),[-ampl,-ampl],color='b',ls ='--')
#plt.annotate(r'$A = $'+repr(ampl),fontsize = 20, xy=(25, -ampl), xytext=(2, -2.5), arrowprops=dict(facecolor='b', shrink=0.05),)
#plt.annotate(r'$A_{MIA} = $'+repr(ampl)+r'$\cdot \sqrt{\frac{\beta1_{MIA}}{\beta2_{MIA} }}$',fontsize = 20, xy=(37, -ampl*coefm), xytext=(27, -2.3), arrowprops=dict(facecolor='r', shrink=0.05),)
text1 = plt.text(5,2.3,r'$\nu_x = 18.84$'+'\n'+r'$\nu_y = 18.73$',fontsize = 18)

plt.grid(True)
axes = plt.gca()
axes.set_ylim([-3.5,3.5])
plt.legend(loc='best',frameon=True)

plt.title(r'Beam history on BPMs')
plt.xlabel(r'n, turn', fontsize = 18)
plt.ylabel(r'x, mm',fontsize = 18)

#save(name='track_first_wp')
'''
#------------------------PLOT--------------NUMBER-OF-TURNS-----------------------------------------------------------------
'''
fig = plt.figure()

plt.subplot(211)
gr2 = plt.plot(per,betaerror, color = 'b', label = r'Error of $\beta$ calculations')
gr4 = plt.plot(per,phaserror, color = 'r',label = r'Error of $\Psi$ calculations')

text1 = plt.text(50,14,r'$\nu_x = 18.84$'+'\n'+r'$\nu_y = 18.73$',fontsize = 18)

plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,241])
axes.set_ylim([0,20])
plt.legend(loc='best',frameon=True, fontsize = 15)

plt.title(r'MIA error dependence on the length of beam history',fontsize = 17)
plt.ylabel(r'$\delta$, %',fontsize = 18)

plt.subplot(212)

gr2 = plt.plot(per,betaerror2, color = 'b', label = r'Error of $\beta$ calculations')
gr4 = plt.plot(per,phaserror2, color = 'r',label = r'Error of $\Psi$ calculations')
text1 = plt.text(50,14,r'$\nu_x = 17.42$'+'\n'+r'$\nu_y = 17.37$',fontsize = 18)
text1 = plt.text(120,17,r'noise = 100 $\mu m$')
#plt.annotate(r'$\sim \frac{|cos(2\cdot 2 \pi \nu_x n)|}{n}$',fontsize = 25, xy=(100, 5), xytext=(140, 8), arrowprops=dict(facecolor='black', shrink=0.05),)

plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,241])
axes.set_ylim([0,20])
#plt.legend(loc='best',frameon=True)

plt.xlabel(r'N, length of tracks')
plt.ylabel(r'$\delta$, %')

#save(name='te')
'''
#----------------------------------PLOT---------NOISE---------------------------------------------------
'''
per = np.array(range(0, 300, 5))
fig = plt.figure()

plt.subplot(211)
gr1 = plt.plot(0.001*per/A,betaerror, color = 'r', label = r'Error of $\beta$, A = 2 $mm$')
gr2 = plt.plot(0.001*per/A,phaserror, color = 'b',label = r'Error of $\Psi$, A = 2 $mm$')
gr3 = plt.plot(0.001*per/A,betaerror12,ls = '--', color = 'r')
gr4 = plt.plot(0.001*per/A,phaserror12, ls = '--',color = 'b')

text1 = plt.text(0.11,2.5,r'200 turns')
text1 = plt.text(0.08,2,r'$\nu_x = 18.84$'+'\n'+r'$\nu_y = 18.73$')

plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,0.146])
axes.set_ylim([0,3])
plt.legend(loc='upper left',frameon=True, fontsize = 15)

plt.title(r'MIA error dependence on noise')
plt.ylabel(r'$\delta$, %')

plt.subplot(212)
gr1 = plt.plot(0.001*per/A,betaerror2, color = 'r')
gr2 = plt.plot(0.001*per/A,phaserror2, color = 'b')
gr3 = plt.plot(0.001*per/A,betaerror22,ls = '--', color = 'r', label = r'Error of $\beta$,  A = 4 $mm$')
gr4 = plt.plot(0.001*per/A,phaserror22,ls = '--', color = 'b',label = r'Error of $\Psi$, A = 4 $mm$')

text1 = plt.text(0.08,2,r'$\nu_x = 17.42$'+'\n'+r'$\nu_y = 17.37$')

plt.grid(True)
axes = plt.gca()
axes.set_xlim([0,0.146])
axes.set_ylim([0,3])
plt.legend(loc='upper left',frameon=True,fontsize = 15)

plt.xlabel(r'$noise$, $\Delta x/A$')
plt.ylabel(r'$\delta$, %')

#save(name='ne')
'''

plt.show()
#----------------------------------------------------------------------------------------------------------


f.close()
