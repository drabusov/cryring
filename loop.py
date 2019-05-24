import os
import numpy as np
import matplotlib.pyplot as plt

def calculate_k(q1,q2):
	kf= 0.5086546699 + (q1-2.42)*0.192535 + (q2-2.42)*0.0939387
	kd= -0.6511149282 -(q1-2.42)*0.178709 -(q2-2.42)*0.34252
	return kf,kd	

def change_parameter(p1,p2):

	f =open("Optics_v7.str")
	g = open("optics.str", "w")


	[g.write(line) for line in f if "kq" not in line]

	g.write("\nkqfl := {};".format(p1))
	g.write("\nkqdl := {};".format(p2))	

	f.close()
	g.close()


def read_tunes():
	
	f = open("twiss.txt")
	for line in f:
		line = line.strip()
		line = line.split()
		if line[1] == "Q1":
			Q1 = float(line[3])
		if line[1] == "Q2":
			Q2 = float(line[3])
			break

	return Q1,Q2


def read_beta():
    
	twiss = open("twiss.txt")
	header = 48

	s, sbpm, beta, betabpm = list(), list(), list(), list()

	for i,line in enumerate(twiss):
		line=line.strip()
		line=line.split()
		if line[1] == "Q1":
			Q1 = float(line[3])
		if line[1] == "Q2":
			Q2 = float(line[3])
		if i > header:
			s.append(float(line[2]))
			beta.append(float(line[6]))
			if line[1]=='"HMONITOR"':
				sbpm.append(float(line[2]))
				betabpm.append(float(line[6]))


	twiss.close()
	return Q1, Q2, s, sbpm, beta, betabpm

read_beta()

#l_arr = [(0.4884152, -0.6348821), (0.4864899,-0.633095), (0.484565, -0.631308), (0.4826392, -0.6295209), (0.4807138, -0.6277338)]
#l_arr = [(x, -x-0.142460312) for x in np.linspace(0.4807138, 0.5186546699,200)]
l_arr = [(2.3,2.43), (2.32,2.45),(2.32,2.46),(2.32,2.47)]

for l in l_arr:
	a =calculate_k(*l)
	print(a)

'''
f = open("data.txt", "w")
#f.write("Q1, Q2, beta_mean, beta_bpm_mean\n")	
f.write("kqfl, kqdl, Q1, Q2\n")	


for l in l_arr:
	change_parameter(*l)
	os.system("./madx cryring_ideal.madx -> log.txt")
#	Q1, Q2, s, sbpm, beta, betabpm =read_beta()
#	level1 = np.mean(beta)
#	level2 = np.mean(betabpm)
#	f.write("{}, {}, {}, {}\n".format(Q1,Q2,level1,level2))	

	Q1,Q2 =read_tunes()

	f.write("{}, {}, {}, {}\n".format(*l, Q1,Q2))	
'''


'''
l_arr = np.linspace(0.285,0.3, 100)


q1, q2 = [], []

for l in l_arr:
	change_parameter(l)
	os.system("./madx cryring_track.madx -> log.txt")
	a,b =read_tunes()
	q1.append(a)
	q2.append(b)
#print("{} {}".format(a,b))

plt.figure()
plt.plot(l_arr,q1, label = "horizontal tune")
plt.plot(l_arr,q2, label = "vertical tune")
plt.xlim(min(l_arr),max(l_arr))
plt.title("Tune vs length of quad")
plt.xlabel("Length, m")
plt.ylabel("tune")
plt.legend()
plt.grid(True)
plt.savefig("tune_vs_len.pdf")

plt.show()
'''
