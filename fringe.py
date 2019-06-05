import os
import numpy as np
from scipy.optimize import least_squares as ls

import matplotlib.pyplot as plt

#----------------------------------------------------------------------------

def read_twiss(filename, header):
    
    twiss = open(filename)  # This file updates every time of madx running

    betaxbpm, betaybpm = list(), list()

    for i,line in enumerate(twiss):
        line=line.strip()
        line=line.split()
        if i > header:
            if line[0] in bpms_x:
                betaxbpm.append(float(line[2]))

            if line[0] in bpms_y:
                betaybpm.append(float(line[3]))

    twiss.close()
    return betaxbpm, betaybpm

#----------------------------------------------------------------------------

#corr = ['mhq','mhs', 'mqs', 'mqo']
corr = ['mhs']
theta0 = [np.random.normal(0,0.0001) for elem in corr] # vector of parameters to optimize, first point to start

# 29.03
measured_x =[0.0003309073729546661, 0.0002212738428333959, 0.00010740585805192146, 0.0002516754207025145, 0.0002681394558083527, 6.721085594107403e-05, 0.00032595340641084424]
measured_y = [7.873070365921625e-05, 0.00010167789015970228, 3.648461631362319e-05, 8.931412173066063e-05, 7.680162282391101e-05, 3.5642186356593414e-05, 7.21851181894217e-05]

# 8.04
#measured_x = [0.0009299879954459251, 0.0008033671774755714, 0.0002483917301924922, 0.0007132292204245833, 0.0006288374903291293, 0.00026724949197784405, 0.0007892775300489714]
#measured_y = [0.0001403398105088781, 0.00020603624848912746, 7.13593115559004e-05, 0.00019439871384627867, 0.00017004167603238752, 5.205934870703545e-05, 0.0001329893543921361]


bpms_x = ['"YR02DX1"', '"YR06DX1"', '"YR07DX1"','"YR08DX1"','"YR10DX1"','"YR11DX1"','"YR12DX1"']
bpms_y = ['"YR02DX2"', '"YR06DX2"','"YR07DX2"','"YR08DX2"','"YR10DX2"','"YR11DX2"','"YR12DX2"']


def metric_counter(theta):

    f = open('Optics_orbit_first.str')
#    f = open('Optics_orbit_second.str')
    g = open('Optics_test.str', 'w')

    for line in f:
        if line!='\n':
            line = line.rstrip(';\n')
            line = line.split(':=')
            if not line[0].strip('\t') in corr:
                g.write('{}:={};\n'.format(line[0],line[1]))


    for i,v in enumerate(theta):
        g.write('{}\t:={};\n'.format(corr[i], v))


    f.close()
    g.close()


    os.system("./madx cryring_fringe.madx -> out.dat")

    SP2 = 48 # header of twiss file
    betaxbpm, betaybpm = read_twiss('twiss.txt',SP2)



    metric_x = (np.array(measured_x/np.mean(measured_x)) - np.array(betaxbpm/np.mean(betaxbpm)))
    metric_y = (np.array(measured_y/np.mean(measured_y)) - np.array(betaybpm/np.mean(betaybpm)))
#    metric_orb = [np.std(orbit_x),np.std(orbit_y)]



    print("metric is {}".format(np.mean(metric_x)+np.mean(metric_y)))

    return list(metric_x)+list(metric_y)
#    return list(metric_x)+list(metric_y)+metric_orb




ls(metric_counter, theta0)

#metric_counter(theta0)


SP2 = 48 # header of twiss file
betaxbpm, betaybpm = read_twiss('twiss.txt',SP2)

plt.figure()
plt.plot(betaxbpm/np.mean(betaxbpm), color = "red", ls="--")
plt.plot(measured_x/np.mean(measured_x), marker = "^", color = "red")

plt.plot(betaybpm/np.mean(betaybpm), color = "blue", ls="--")
plt.plot(measured_y/np.mean(measured_y), marker = "^", color = "blue")

plt.grid(True)


plt.show()









