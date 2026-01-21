
#* Mathematical imports
import math
import numpy as np

#* Utility imports
import os
import sys
import csv
import shutil
import datetime
import logging
import subprocess
import configparser
from subprocess import call
from decimal import Decimal

#* Matplotlib imports
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib.collections import LineCollection
from matplotlib.ticker import FormatStrFormatter

#* Addittional options for imported libraries
mpl.use('Agg')
mpl.rcParams['agg.path.chunksize'] = 100000
plt.rc('axes', labelsize=18)
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

#* Additional Functions
#* Locate CBwaves binary
def find(name, path):
	for root, dirs, files in os.walk(path):
		if name in files:
			return os.path.join(root, name)

#! If you want to remove Data files after simulations ran and figire was made (see below).
removedat = False

#! If you want to run your simulation after the Last stable orbit
after6M = True

#* Definition of auxiliary variables

pi = 3.141592653589793		#* Pi
c = 2.99792458e8			#* Speed of light in [m/s]
G = 6.674184e-11			#* Gravitational constant 
Gc2 = G/(c*c)				#* G/c^2
ly = 9460730472580800.0		#* Light year
pc = 3.26156 * ly			#* Parsec
year = 365.25636*24.*3600.	#* Orbital period time of Earth in seconds
msun = 1476.62504			#* Mass of the sun in G=c=1 system, [m]

#! - Where is the cbwaves executable
if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
	path = os.path.dirname(os.getcwd())
	cbwbinary = find('cbwaves',path)
elif sys.platform.startswith('win32'):
	path = os.path.dirname(os.getcwd())
	cbwbinary = find('cbwaves',path)

#! - Output file name prefix
filenameprefix = "OJ287"

#! - Mass ratio
M1 = 1.8348*math.pow(10, 10)
M2 = 1.5013*math.pow(10, 8)
m1 = M1*msun
m2 = M2*msun
M = m1 + m2

#! - orbital time
T = 12.062 * year
Tc = T*c

#! - orbit frequency
f = 1 / T

#! - Initial separation:
z = 0.3
r0 = 16.41369744591233
r = r0 * (M * (1 + z))
T = 2. * pi * r/(c * math.sqrt(M/r))

#! - Time step
dt = 100000

#! - Maximum evolution time
tmax = 3.0e20

#! - Maximum number of orbits
orbitsmax = 1e6

#! - Eccentricity of the orbit
epsilon = 0.657

#! - Polar angles
iota = 0.0
phi = 0.0
theta = 0.0
varphi = 0.0
psi = 0.0

#! - Spin definition
s1 = 0.381
s2 = 0
s1x = 0
s1y = 0
s1z = s1*math.cos(0)
s2x = 0
s2y = 0
s2z = 0

#! - Output filename definition
outfile = "cbwaves.out"
ftfile = "cbwaves.ft"

#! - Active corrections for the center of mass acceleration terms and for the radiation field
corrs = "'PN','2PN','SO','SS','RR','PNSO','3PN','1RR','2PNSO','RRSO','RRSS','4PN'"
hterms = "'Q','P05Q','PQ','P15Q','P15Qtail','PQSO','P15QSO','P2Q','PQSS'"

#! - Output variables
outvars = "t,x1,y1,z1,x2,y2,z2,h,s1x,s1y,s1z,E_tot,E_rad,h_+,h_x,orbfreq,r,ecc_r,v2,mr,orbits"

#! All the possible output variables:
#* t,rNx,rNy,rNz,vNx,vNy,vNz,rPNx,rPNy,rPNz,vPNx,vPNy,v2PNz,r2PNx,r2PNy,r2PNz,v2PNx,v2PNy,v2PNz,
#* r3PNx,r3PNy,r3PNz,v3PNx,v3PNy,v3PNz,r4PNx,r4PNy,r4PNz,v4PNx,v4PNy,v4PNz,rSOx,rSOy,rSOz,vSOx,
#* vSSy,vSSz,rSSx,rSSy,rSSz,vSSx,vSSy,vSSz,rRRx,rRRy,rRRz,vRRx,vRRy,vRRz,
#* rPNSOx,rPNSOy,rPNSOz,vPNSOx,vPNSOy,vPNSOz,r1RRx,r1RRy,r1RRz,v1RRx,v1RRy,v1RRz,
#* rRRSOx,rRRSOy,rRRSOz,vRRSOx,vRRSOy,vRRSOz,rRRSSx,rRRSSy,rRRSSz,vRRSSx,vRRSSy,vRRSSz,
#* s1x,s1y,s1z,s2x,s2y,s2z,orbits,rx,ry,rz,r,mr,vx,vy,vz,v,v2,h_+,h_x,h,E,E_N,E_PN,E_2PN,E_3PN,
#* E_4PN,E_SO,E_SS,E_PNSO,E_rad,E_PNrad,E_2PNrad,E_SOrad,E_SSrad,E_Nrad,E_25PNrad,E_PNSOrad
#* E_PNtot,E_2PNtot,E_SOtot,E_SStot,E_tot,Jx_Nrad,Jy_Nrad,Jz_Nrad,J_tot,Jx,Jy,Jz,J2,x1,y1,z1
#* x2,y2,z2,hkp,hp22,hx22,hp21,hx21,hp2m2,hx2m2,hp2m1,hx2m1,hp20,hx20,ecc_N,ecc_PN, ecc_ell,ecc_r
#* ecc_RL,orbits_turn,Lnx,Lny,Lnz,Ln2,Lx,Ly,Lz,L2,,drAdiab,rdot,orbfreq
#! for refernce see cbwaves.cxx line 1657 to 1819

#! - Do we want initial eccentricity approximation
eccapprox = "no"

#! - Do we want checkpointing
checkpoint = "no"

#! - What is the description of the run - to be used for the checkpoint files
description = "OJ287 orbital and waveform analysis in 4PN"

#! - distance from observer: D = 2*mu = 2*m1*m2/(m1+m2)
D = 1.64689e9*pc

#! - logging level
loglevel = 6

#! - simulations stops when r<rmin or r>rmax
rmin = 3.*(m1 + m2)
rmax = 150*(m1 + m2)

#! - loop variables
i = []
j = []
k = []

#! print variables
printstep = 1
printorbit = 0
adaptive = "no"
adaptive_step = 1000


#! - gauge parameters for RR
eta = m1 * m2 / (m1 + m2) / (m1 + m2)
galpha = 4
gbeta = 5
gdelta1 = -99/14.+27.*eta
gdelta2 = 5.*(1. - 4.*eta)
gdelta3 = 274./7. + 67./21.*eta
gdelta4 = 5./2.*(1. - eta)
gdelta5 = -1./7.*(292. + 57.*eta)
gdelta6 = 51./28. + 71./14.*eta

#!
#! - For loops comes here
#!

inifile = filenameprefix + '.ini'
outfile = filenameprefix + '.dat'

logging.info("OJ 287 simulation is runing.")

##!
#! - Creating the .ini file
##!
if not os.path.isfile(inifile):
	config = configparser.RawConfigParser()
	config.optionxform = str
	config['output'] = {

		'outfile': outfile,
		'ftfile': ftfile,
		'outvars': outvars,

	}
	config['input'] = {

		'm1': m1,	                             # mass of bigger star
		'm2': m2,	                             # mass of lighter star
		'tmax': tmax,                            # maximum evolution time
		'orbitsmax': orbitsmax,                  # maximum number of revolutions
		'T': T,		                             # orbit time
		'f': f,                                  # orbiting frequency
		'dt': dt,                                # RK4 time step
		'epsilon': epsilon,                      # eccentricity
		'rmin': rmin,                            # final relative distance
		'rmax': rmax,                            # final relative distance
		'r': r,                                  # initial separation
		# 'ri': ri,                                # initial separation for hyperbolic orbits
		'D': D,                                  # distance to observer
		'iota': iota,                            # polar angle in the source frame
		'phi': phi,                              # polar angle in the source fram
		'theta': theta,                          # Euler angle between TT and detector frame
		'varphi': varphi,                        # Euler angle between TT and detector frame
		'psi': psi,                              # Euler angle between TT and detector frame
		's1x': s1x,                              # spin x components of m1
		's1y': s1y,                              # spin y components of m1
		's1z': s1z,                              # spin z components of m1
		's2x': s2x,                              # spin x components of m2
		's2y': s2y,                              # spin y components of m2
		's2z': s2z,                              # spin z components of m2
		'hterms': hterms,                        # PN order for the waveform
		'corrs': corrs,                          # PN order for the motion
		# approximate or not the initial eccentricity
		'eccapprox': eccapprox,
		'checkpoint': checkpoint,                # wheter to make or not checkpoint files
		'description': description,              # human readable description of the run
		'printstep': printstep,                  # print variables after every X step
		'printorbit': printorbit,                # print variables after every X orbit
		'loglevel': loglevel,                    # logging level (0 - 6)
		'adaptive': adaptive,                    # apply adaptive time step or not
		'adaptive_step': adaptive_step,          # how many time step per revolution
		'alpha': galpha,                         # newtonian radiation term gauge parameter
		'beta': gbeta,                           # newtonian radiation term gauge parameter
		# post-newtonian radiation term gauge parameter
		'delta1': gdelta1,
		# post-newtonian radiation term gauge parameter
		'delta2': gdelta2,
		# post-newtonian radiation term gauge parameter
		'delta3': gdelta3,
		# post-newtonian radiation term gauge parameter
		'delta4': gdelta4,
		# post-newtonian radiation term gauge parameter
		'delta5': gdelta5,
		# post-newtonian radiation term gauge parameter
		'delta6': gdelta6,
	}
	with open(inifile, 'w') as configfile:
		config.write(configfile)
else:
	logging.info("INI file already exist. Starting next process!")

#! Gives the configuration files (e.g test.ini) to CBwaves binary
if not os.path.isfile(outfile):
	call([cbwbinary, inifile])
	logging.info("The data file is generated: %s", outfile)
else:
	logging.info("Data file already exist. Starting next process!")

#!
#! Some example plots is created her
#!

logging.info("Starts to create example plots!")

#! Pre calculations for the plots
A = (G*M*eta)/D

###

tt = []			# row[0]
xx1 = []		# row[1]
yy1 = []		# row[2]
zz1 = []		# row[3]
xx2 = []		# row[4]
yy2 = []		# row[5]
zz2 = []		# row[6]
hh = []			# row[7]
s1x = []		# row[8]
s1y = []		# row[9]
s1z = []		# row[10]
etot = []		# row[11]
erad = []		# row[12]
hhp = []		# row[13]
hhx = []		# row[14]
orbfreq = []	# row[15]
rr = []			# row[16]
eccr = []		# row[17]
vv2 = []		# row[18]
mmr = []		# row[19]


with open(outfile) as zw:
	zwval = csv.reader(zw, delimiter=' ')
	for row in zwval:
		tt.append(float(row[0]))
		xx2.append(float(row[4]))
		yy2.append(float(row[5]))
		zz2.append(float(row[6]))

rwb = 0
# yearh = 365.25636*24.
x2 = np.array([i/M for i in xx2])
y2 = np.array([j/M for j in yy2])
z2 = np.array([k/M for k in zz2])
# th = np.array([i/3600 for i in tt])
ty = np.array([i/year for i in tt])

#!
#! Colored plot of m2's orbits
#!
logging.info("Ploting after LSO orbits of the auxilary mass \$m2\$.")
if after6M:
	x2_6m = x2[5482000:]
	y2_6m = y2[5482000:]
	z2_6m = z2[5482000:]
	ty_6m = ty[5482000:]

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	norm = plt.Normalize(ty_6m.min(), ty_6m.max())
	points = np.array([x2_6m, y2_6m, z2_6m]).T.reshape(-1, 1, 3)
	segments = np.concatenate([points[:-1], points[1:]], axis=1)

	lc3d = Line3DCollection(segments, linewidths=1, cmap='jet')
	lc3d.set_array(ty_6m)
	line = ax.add_collection(lc3d)
	fig.colorbar(line, ax=ax, label='t [year]')

	ax.set_xlim(x2_6m.min(), x2_6m.max())
	ax.set_ylim(y2_6m.min(), y2_6m.max())
	ax.set_zlim(z2_6m.min(), z2_6m.max())

	ax.set_xlabel('x/M')
	ax.set_ylabel('y/M')
	ax.set_zlabel('z/M')
	ax.grid(visible=True, which='major', axis='both')

	fig.tight_layout()
	plt.savefig('ColoredOrbits6M_' + filenameprefix + '.png')
	plt.close()

logging.info("Ploting orbits of the auxilary mass \$m2\$.")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

norm = plt.Normalize(ty.min(), ty.max())
points = np.array([x2, y2, z2]).T.reshape(-1, 1, 3)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

lc3d = Line3DCollection(segments, linewidths=1, cmap='jet')
lc3d.set_array(ty)
line = ax.add_collection(lc3d)
fig.colorbar(line, ax=ax, label='t [year]')

ax.set_xlim(x2.min(), x2.max())
ax.set_ylim(y2.min(), y2.max())
ax.set_zlim(z2.min(), z2.max())

ax.set_xlabel('x/M')
ax.set_ylabel('y/M')
ax.set_zlabel('z/M')
ax.grid(visible=True, which='major', axis='both')

fig.tight_layout()
plt.savefig('ColoredOrbits_' + filenameprefix + '.png')
plt.close()

del x2, y2, z2, xx1, yy2, zz2

#!
#! Energy plots
#!
logging.info("Ploting the total and radiated energies of the system")
with open(outfile) as zw:
	zwval = csv.reader(zw, delimiter=' ')
	for row in zwval:
		etot.append(float(row[11]))
		erad.append(float(row[12]))
		rr.append(float(row[16]))

r = np.array([i/M for i in rr])
Etot = np.array([i for i in etot])
Erad = np.array([i for i in erad])

if after6M:
	r6M = r[5482000:]
	Etot6M = Etot[5482000:]
	Erad6M = Erad[5482000:]
	ty_6m = ty[5482000:]

	fig = plt.figure('5b')
	fig.clf()
	ax1 = fig.add_subplot(2, 1, 1)
	cmap1 = cm.jet
	sc1 = ax1.scatter(ty_6m, Etot6M, label='$E_{tot}$', c=r6M, s=0.1, cmap=cmap1)
	plt.xlabel('t [year]')
	plt.ylabel('E')
	ax1.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
	plt.legend()
	plt.grid(visible=True, which='both')
	cbar = plt.colorbar(sc1)
	cbar.set_label('r/M')

	ax2 = fig.add_subplot(2, 1, 2)
	cmap2 = cm.gnuplot
	sc2 = ax2.scatter(ty_6m, Erad6M, label='$E_{rad}$', c=r6M, s=0.1, cmap=cmap2)
	plt.xlabel('t [year]')
	plt.ylabel('E')
	ax2.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
	plt.legend()
	plt.grid(visible=True, which='both')
	cbar2 = plt.colorbar(sc2)
	cbar2.set_label('r/M')

	fig.tight_layout()
	plt.savefig('Energy6M_' + filenameprefix + '.png')
	plt.close()


fig = plt.figure('5b')
fig.clf()
ax1 = fig.add_subplot(2, 1, 1)
cmap1 = cm.hsv
sc1 = ax1.scatter(r, Etot, label='$E_{tot}$', c=ty, s=0.1, cmap=cmap1)
plt.xlabel('r/M')
plt.ylabel('E')
ax1.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
ax1.set_xlim(100,0)
plt.legend()
plt.grid(visible=True, which='both')
cbar = plt.colorbar(sc1)
cbar.set_label('t [year]')

ax2 = fig.add_subplot(2, 1, 2)
cmap2 = cm.twilight_shifted
sc2 = ax2.scatter(r, Erad, label='$E_{rad}$', c=ty, s=0.1, cmap=cmap2)
plt.xlabel('r/M')
plt.ylabel('E')
ax2.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
ax2.set_xlim(100,0)
plt.legend()
plt.grid(visible=True, which='both')
cbar2 = plt.colorbar(sc2)
cbar2.set_label('t [year]')

fig.tight_layout()
plt.savefig('Energy_' + filenameprefix + '.png')
plt.close()

del etot, erad, Etot, Erad

#!
#! Waveform plots
#!
logging.info("Ploting the Gravitational Waves emitted by the system.")
with open(outfile) as zw:
	zwval = csv.reader(zw, delimiter=' ')
	for row in zwval:
		hh.append(float(row[7]))
		hhp.append(float(row[13]))
		hhx.append(float(row[14]))

h = [i for i in hh]
hp = [j for j in hhp]
hx = [k for k in hhx]

if after6M:
	h6M = h[5482000:]
	hp6M = hp[5482000:]
	hx6M = hx[5482000:]
	r6M = r[5482000:]
	ty_6m = ty[5482000:]

	fig = plt.figure('Waveform')
	plt.plot(ty_6m, h6M, 'r', label='Waveform')
	plt.xlabel('t [year]')
	plt.ylabel('$h$')
	plt.grid(visible=True, which='both')
	fig.tight_layout()
	plt.savefig('Waveform(6M)_' + filenameprefix + '.png')
	plt.close()

	fig = plt.figure('Waveform(+)')
	plt.plot(ty_6m, hp6M, 'r', label='Waveform')
	plt.xlabel('t [year]')
	plt.ylabel('${h^{\\theta\\theta}}_{TT}$')
	plt.grid(visible=True, which='both')
	fig.tight_layout()
	plt.savefig('Waveform(+,6M)_' + filenameprefix + '.png')
	plt.close()

	fig = plt.figure('Waveform(x)')
	plt.plot(ty_6m, hx6M, 'r', label='Waveform')
	plt.xlabel('t [year]')
	plt.ylabel('${h^{ \\theta \phi}}_{TT}$')
	plt.grid(visible=True, which='both')
	fig.tight_layout()
	plt.savefig('Waveform(x,6M)_' + filenameprefix + '.png')
	plt.close()

fig = plt.figure('Waveform')
plt.plot(ty, h, 'r', label='Waveform')
plt.xlabel('t [year]')
plt.ylabel('$h$')
plt.grid(visible=True, which='both')
fig.tight_layout()
plt.savefig('Waveform_' + filenameprefix + '.png')
plt.close()

fig = plt.figure('Waveform(+)')
plt.plot(ty, hp, 'r', label='Waveform')
plt.xlabel('t [year]')
plt.ylabel('${h^{\\theta\\theta}}_{TT}$')
plt.grid(visible=True, which='both')
fig.tight_layout()
plt.savefig('Waveform(+)_' + filenameprefix + '.png')
plt.close()

fig = plt.figure('Waveform(x)')
plt.plot(ty, hx, 'r', label='Waveform')
plt.xlabel('t [year]')
plt.ylabel('${h^{ \\theta \phi}}_{TT}$')
plt.grid(visible=True, which='both')
fig.tight_layout()
plt.savefig('Waveform(x)_' + filenameprefix + '.png')
plt.close()

del h, hp, hx, hh, hhp, hhx

#!
#! Spin's components
#!
logging.info("Ploting the components of the spin.")
with open(outfile) as zw:
	zwval = csv.reader(zw, delimiter=' ')
	for row in zwval:
		s1x.append(float(row[8]))
		s1y.append(float(row[9]))
		s1z.append(float(row[10]))

fig = plt.figure('s1z')
plt.plot(ty, s1z, label='Waveform')
plt.xlabel('t [year]')
plt.ylabel('$S^{(1)}_z$')
plt.grid(visible=True, which='both')
plt.savefig('Spin1_z.png')
plt.close()

del s1x, s1y, s1z

#!
#! Plot of eccentricity as function of t
#!
logging.info("Ploting thr eccentricity as fuction of t.")
with open(outfile) as zw:
	zwval = csv.reader(zw, delimiter=' ')
	for row in zwval:
		eccr.append(float(row[17]))

ecc = np.array([i for i in eccr])

fig = plt.figure()
ax = fig.gca()

norm = plt.Normalize(r.min(), r.max())
points = np.array([ty, ecc]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

lc = LineCollection(segments=segments, linewidths=2, cmap='jet')
lc.set_array(r)
line = ax.add_collection(lc)
fig.colorbar(line, ax=ax, label='r/M', pad=0)

ax.set_xlim(ty.min(), ty.max())
ax.set_ylim(0, 1.1)
ax.set_xlabel('t[year]')
ax.set_ylabel('excentrictás')

ax.grid(visible=True, which='major', axis='both')
fig.tight_layout()
plt.savefig('Eccentricity_'+ filenameprefix +'.png')
plt.close()

del ecc, eccr


#! cleaning up generated files
now = datetime.datetime.now()
dirname = 'data/' + filenameprefix + '_' + now.strftime("%y-%m-%d-%H%M")

if not os.path.isdir('data'):
	call(['mkdir','data'])

call(['mkdir', dirname])
logging.info(dirname + " has been created!")

if os.path.isdir(dirname):
    os.system('mv *.png ' + dirname)
    os.system('mv *.ini ' + dirname)

    if removedat:
        os.system('rm *.dat')
        logging.info("Generated data files has been removed!")
    else:
       os.system('mv *.dat ' + dirname)
       logging.info('Generated data file has been moved to' + dirname + '!')
else:
    print("Directory doesn\'t exist!")
