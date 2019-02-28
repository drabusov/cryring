import os
import numpy as np
from scipy.optimize import least_squares as ls

#----------------------------------------------------------------------------

def read_twiss(filename, header):
    
    twiss = open(filename)  # This file updates every time of madx running

    s, sbpm, beta, betabpm, phasebpm = list(), list(), list(), list(), list()

    for i,line in enumerate(twiss):
        line=line.strip()
        line=line.split()
        if i > header:
            s.append(float(line[2]))
            beta.append(float(line[6])) 
            if 'DX1' in line[0]:
                sbpm.append(float(line[2]))
                betabpm.append(float(line[6]))
                phasebpm.append(2*np.pi*float(line[10]))
    twiss.close()
    return betabpm

#----------------------------------------------------------------------------

corr = ['qf1','qf3','qf5','qd2','qd4','qd6']
theta0 = [np.random.normal(0,0.01) for elem in corr] # vector of parameters to optimize, first point to start

def metric_counter(theta):

    f = open('Optics_v7.str')
    g = open('Optics_test.str', 'w')

    for line in f:
        if line!='\n':
            line = line.rstrip(';\n')
            line = line.split(':=')
            if 'kqfl' in line[0]:
                kqfl = float(line[1])
            if 'kqdl' in line[0]:
                kqdl = float(line[1])
            g.write('{} := {};\n'.format(line[0],line[1]))

    for i,v in enumerate(theta):
        if 'qf' in corr[i]:
            g.write('{} := {};\n'.format(corr[i],kqfl+v))
        else:
            g.write('{} := {};\n'.format(corr[i],kqdl+v)) 

    f.close()
    g.close()

    os.system("/home/laptop/cryring/madx < cryring.madx > ./out.dat")


    SP2 = 48 # header of twiss file
    betabpm = read_twiss('twiss.txt',SP2)

    index = [0,2,4,5,7]
    beta = [x for i,x in enumerate(betabpm) if i in index]
    metric = np.std(beta)/np.mean(beta)

    print('periodic beta is {}'.format(beta))
    print('set of parameters is {}'.format(theta))
    print('metric is {}'.format(metric))


    return metric**2

ls(metric_counter, theta0, xtol = 10**(-5))

