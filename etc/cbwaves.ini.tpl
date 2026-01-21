#
# Initialization file template for the CBwaves executable
# 

[output]
outfile = /tmp/mywaveform.dat            # output filename
ftfile = /tmp/mywaveform.ft              # file for the Fourier-transformed
outvars = "t,orbits,rx,ry,rz,h_+,E_N,E_PNtot,E_tot,E_rad,hp22,hx22,x1,y1,z1,x2,y2,z2,hp22,hp21,hp20,hp2m1,hp2m2,hx22,hx21,hx20,hx2m1,hx2m2,h_x,h"
                                         # quantitites to output

[input]
m1        = 26730                        # mass of first star in meters
m2        = 14850                        # mass of lighter star in meters
tmax      = 10e5                         # maximum evolution time in sec
orbitsmax = 10e5                         # maximum number of revolutions
T         = 0.05                         # initial orbit time
f         = 18                           # initial orbiting frequency
dt        = 1e-05                        # initial RK4 time step
epsilon   = 0                            # initial eccentricity
rmin      = 101574.6                     # final alloved mimimum separation
rmax      = 9999999999999                # final allowed maximum separation (for open orbits)
r         = 304041.0                     # initial separation
D         = 2.465e+23                    # distance to observer
iota      = 0                            # inclination angle of orbit
phi       = 0                            # polar angle in source frame
theta     = 0                            # polar angle in source frame
varphi    = 0                            # -- angle
psi       = 0                            # -- angle
s1x       = 0                            # spin1 x components of m1
s1y       = 0                            # spin1 y components of m1
s1z       = 0                            # spin1 z components of m1
s2x       = 0                            # spin2 x components of m2
s2y       = 0                            # spin2 y components of m2
s2z       = 0                            # spin2 z components of m2
hterms    = 'Q','P05Q','PQ','P15Q','P15Qtail','PQSO','P15QSO','P2Q','PQSS'
                                         # PN order of calculation for the waveform
corrs     = 'PN','2PN','SO','SS','RR','PNSO','3PN','1RR','2PNSO','RRSO','RRSS'
                                         # PN order for the movement of the bodies
printstep = 1                            # print variables after every X hintegration step
printorbit = 1                           # print variables after every X orbit
checkpoint= yes                          # do you want checkpoint or not
description=testrun                      # description of the run
loglevel  = 4                            # level of verbosity
adaptive   = no
adaptive_step = 100
alpha  = 4;                              # newtonian radiation term gauge parameters
beta   = 5;
delta1 = -0.872448979591837;             # post-newtonian radiation term gauge parameters
delta2 = 0.408163265306123;
delta3 = 39.8753644314869;
delta4 = 1.92602040816327;
delta5 = -43.5838192419825;
delta6 = 2.98578717201166;