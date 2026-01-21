
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

#* Import export function from export_oj287_coordinates.py
from export_oj287_coordinates import export_coordinates

#* Additional Functions
#* Locate CBwaves binary
def find(name, path):
	for root, dirs, files in os.walk(path):
		if name in files:
			return os.path.join(root, name)

#! If you want to remove Data files after simulations ran and figure was made (see below).
removedat = False

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
cbwbinary = None
if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
	# Look in the parent directory (project root) recursively
	path = os.path.abspath(os.path.join(os.getcwd(), '..'))
	cbwbinary = find('cbwaves', path)
elif sys.platform.startswith('win32'):
	path = os.getcwd()
	cbwbinary = find('CBwaves', path)

if cbwbinary is None:
	logging.error("Could not find 'cbwaves' executable in %s", os.path.abspath(os.path.join(os.getcwd(), '..')))
	sys.exit(1)

#! - Output file name prefix
filenameprefix = "OJ287_precession"

#! ============================================
#! ПАРАМЕТРИ ЗІ СТАТТІ (ApJ Letters 993:L22, 2025, Figure 1)
#! ============================================

#! - Mass ratio q = 0.1 (як у статті, фіолетова лінія)
M1 = 1.8348*math.pow(10, 10)  # Первинна СМЧД в сонячних масах
M2 = M1 * 0.1                 # Вторинна ЧД: q = M2/M1 = 0.1
m1 = M1*msun                   # в метрах (G=c=1)
m2 = M2*msun
M = m1 + m2                    # Повна маса

#! - Initial separation (як у статті: 16.4 r_g)
r0 = 16.4  # бінарна відстань в одиницях r_g
r = r0 * M  # початкова відстань в метрах

#! - orbital time (обчислюється з початкової відстані)
T = 2. * pi * r/(c * math.sqrt(M/r))
Tc = T*c

#! - orbit frequency
f = 1 / T

#! - Time step (ЗМЕНШЕНО для точності прецесії)
dt = 10000  # Зменшено з 100000 до 10000

#! - Maximum evolution time (як у статті: 20,000M ≈ 96 років)
tmax = 96 * year  # ~96 років як у Figure 1 статті

#! - Maximum number of orbits
orbitsmax = 10

#! - Eccentricity of the orbit (реальний ексцентриситет OJ 287)
epsilon = 0.657

#! - Polar angles
iota = 0.0
phi = 0.0
theta = 0.0
varphi = 0.0
psi = 0.0

#! ============================================
#! ORBITAL INCLINATION - УГОЛ НАХИЛУ ОРБІТИ
#! ============================================
#! Цей параметр повертає орбіту від площини x-y
#! 0 градусів = орбіта в площині x-y (як на скріншоті)
#! 30 градусів = орбіта нахилена на 30° від площини x-y
#! 90 градусів = орбіта перпендикулярна до площини x-y
orbital_inclination = 0.0  # градусів

#! Напрямок повороту:
#! 'x' - поворот навколо осі X (нахил в площині y-z)
#! 'y' - поворот навколо осі Y (нахил в площині x-z)
#! 'z' - поворот навколо осі Z (поворот в площині x-y)
rotation_axis = 'x'  # 'x', 'y', або 'z'

#! - Spin definition (ВИМКНЕНО для чистої релятивістської прецесії)
#! Для вивчення чистої 1PN прецесії вимикаємо спін
s1 = 0.0 # Було 0.381, тепер 0 для чистої прецесії
s2 = 0
s1x = 0
s1y = 0
s1z = 0  # s1*math.cos(0) = 0
s2x = 0
s2y = 0
s2z = 0

#! - Output filename definition
outfile = "cbwaves.out"
ftfile = "cbwaves.ft"

#! ============================================
#! ПОПРАВКИ ДЛЯ ЧИСТОЇ РЕЛЯТИВІСТСЬКОЇ ПРЕЦЕСІЇ
#! ============================================
#! Прибрані ВСІ радіаційні поправки (RR, 1RR, RRSO, RRSS)
#! Прибрані спінові поправки (SO, SS, PNSO, 2PNSO) так як спін = 0
#! Залишені тільки PN поправки для чистої прецесії перигелію

corrs = "'PN','2PN','3PN','4PN'"  # Тільки пост-ньютонівські поправки
hterms = "'Q'"  # Мінімальні члени для хвильової форми (не важливо для прецесії)

#! - Output variables (додані змінні для аналізу прецесії)
outvars = "t,x1,y1,z1,x2,y2,z2,r,ecc_r,orbits,E_tot,vx,vy,vz,v2,orbfreq"

#! - Do we want initial eccentricity approximation
eccapprox = "no"

#! - Do we want checkpointing
checkpoint = "no"

#! - What is the description of the run
description = "OJ287 perihelion precession analysis - pure 1PN effect (no radiation, no spin)"

#! - distance from observer: D = 2*mu = 2*m1*m2/(m1+m2)
D = 1.64689e9*pc

#! - logging level
loglevel = 6

#! - simulations stops when r<rmin or r>rmax
#! Збільшено rmax так як без радіації орбіта не буде звужуватися
rmin = 3.*(m1 + m2)
rmax = 500*(m1 + m2)  # Збільшено для стабільності

#! - loop variables
i = []
j = []
k = []

#! print variables (збільшена частота виводу)
printstep = 1
printorbit = 0
adaptive = "no"
adaptive_step = 1000

#! - gauge parameters for RR (не використовуються, але потрібні для ini файлу)
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
#! - Creating INI file and running simulation
#!

inifile = filenameprefix + '.ini'
outfile = filenameprefix + '.dat'
coordinates_file = filenameprefix + '_coordinates.txt'

logging.info("=" * 60)
logging.info("OJ 287 PERIHELION PRECESSION SIMULATION")
logging.info("=" * 60)
logging.info("Studying pure relativistic precession (1PN effect)")
logging.info("Number of orbits: %d", orbitsmax)
logging.info("Eccentricity: %.3f", epsilon)
logging.info("Spin disabled for pure precession study")
logging.info("Radiation reaction DISABLED")
logging.info("=" * 60)

##!
#! - Creating the .ini file
##!
if os.path.isfile(inifile):
	os.remove(inifile)
	logging.info("Removed old INI file")

if os.path.isfile(outfile):
	os.remove(outfile)
	logging.info("Removed old data file")

config = configparser.RawConfigParser()
config.optionxform = str
config['output'] = {
	'outfile': outfile,
	'ftfile': ftfile,
	'outvars': outvars,
}
config['input'] = {
	'm1': m1,
	'm2': m2,
	'tmax': tmax,
	'orbitsmax': orbitsmax,
	'T': T,
	'f': f,
	'dt': dt,
	'epsilon': epsilon,
	'rmin': rmin,
	'rmax': rmax,
	'r': r,
	'D': D,
	'iota': iota,
	'phi': phi,
	'theta': theta,
	'varphi': varphi,
	'psi': psi,
	's1x': s1x,
	's1y': s1y,
	's1z': s1z,
	's2x': s2x,
	's2y': s2y,
	's2z': s2z,
	'hterms': hterms,
	'corrs': corrs,
	'eccapprox': eccapprox,
	'checkpoint': checkpoint,
	'description': description,
	'printstep': printstep,
	'printorbit': printorbit,
	'loglevel': loglevel,
	'adaptive': adaptive,
	'adaptive_step': adaptive_step,
	'alpha': galpha,
	'beta': gbeta,
	'delta1': gdelta1,
	'delta2': gdelta2,
	'delta3': gdelta3,
	'delta4': gdelta4,
	'delta5': gdelta5,
	'delta6': gdelta6,
}
with open(inifile, 'w') as configfile:
	config.write(configfile)
logging.info("INI file created: %s", inifile)

#! Run CBwaves simulation
logging.info("Starting CBwaves simulation...")
call([cbwbinary, inifile])
logging.info("Simulation completed! Data file: %s", outfile)

#! ============================================
#! EXPORT COORDINATES FROM DAT TO TXT
#! ============================================
logging.info("Exporting coordinates from DAT to TXT...")
export_result = export_coordinates(outfile, coordinates_file, max_points=10000)

if export_result is None:
	logging.error("Failed to export coordinates!")
	sys.exit(1)

logging.info("Coordinates exported to: %s", coordinates_file)

#!
#! ============================================
#! READ COORDINATES FROM TXT FILE
#! ============================================
#!

logging.info("Reading coordinates from TXT file for plotting...")

# Arrays for coordinates
ty = []
x1 = []
y1 = []
z1 = []
x2 = []
y2 = []
z2 = []
x_rel = []
y_rel = []
z_rel = []
r = []
orbits = []

with open(coordinates_file, 'r') as f:
	for line in f:
		# Skip header lines (starting with #)
		if line.startswith('#'):
			continue
		
		# Parse data line
		try:
			parts = line.split()
			if len(parts) >= 12:
				ty.append(float(parts[0]))
				x1.append(float(parts[1]))
				y1.append(float(parts[2]))
				z1.append(float(parts[3]))
				x2.append(float(parts[4]))
				y2.append(float(parts[5]))
				z2.append(float(parts[6]))
				x_rel.append(float(parts[7]))
				y_rel.append(float(parts[8]))
				z_rel.append(float(parts[9]))
				r.append(float(parts[10]))
				orbits.append(float(parts[11]))
		except (ValueError, IndexError):
			continue

# Convert to numpy arrays
ty = np.array(ty)
x1 = np.array(x1)
y1 = np.array(y1)
z1 = np.array(z1)
x2 = np.array(x2)
y2 = np.array(y2)
z2 = np.array(z2)
x_rel = np.array(x_rel)
y_rel = np.array(y_rel)
z_rel = np.array(z_rel)
r = np.array(r)
orbits = np.array(orbits)

logging.info("Read %d data points from TXT file", len(ty))

#! ============================================
#! APPLY ORBITAL INCLINATION ROTATION
#! ============================================
if orbital_inclination != 0.0:
    logging.info("Applying orbital inclination: %.1f degrees around %s-axis", orbital_inclination, rotation_axis)
    
    # Convert to radians
    inclination_rad = np.radians(orbital_inclination)
    cos_i = np.cos(inclination_rad)
    sin_i = np.sin(inclination_rad)
    
    if rotation_axis == 'x':
        # Rotation around X-axis (tilts orbit in y-z plane)
        # [1,    0,      0   ]   [x]
        # [0,  cos_i, -sin_i] * [y]
        # [0,  sin_i,  cos_i]   [z]
        x_rel_new = x_rel.copy()
        y_rel_new = y_rel * cos_i - z_rel * sin_i
        z_rel_new = y_rel * sin_i + z_rel * cos_i
        
    elif rotation_axis == 'y':
        # Rotation around Y-axis (tilts orbit in x-z plane)
        # [ cos_i,  0,  sin_i]   [x]
        # [   0,    1,    0  ] * [y]
        # [-sin_i,  0,  cos_i]   [z]
        x_rel_new = x_rel * cos_i + z_rel * sin_i
        y_rel_new = y_rel.copy()
        z_rel_new = -x_rel * sin_i + z_rel * cos_i
        
    elif rotation_axis == 'z':
        # Rotation around Z-axis (rotates in x-y plane)
        # [cos_i, -sin_i, 0]   [x]
        # [sin_i,  cos_i, 0] * [y]
        # [  0,      0,   1]   [z]
        x_rel_new = x_rel * cos_i - y_rel * sin_i
        y_rel_new = x_rel * sin_i + y_rel * cos_i
        z_rel_new = z_rel.copy()
    else:
        logging.warning("Unknown rotation axis '%s', using 'x'", rotation_axis)
        x_rel_new = x_rel.copy()
        y_rel_new = y_rel * cos_i - z_rel * sin_i
        z_rel_new = y_rel * sin_i + z_rel * cos_i
    
    # Apply rotation
    x_rel = x_rel_new
    y_rel = y_rel_new
    z_rel = z_rel_new
    
    logging.info("Orbital inclination applied successfully")

#! ============================================
#! ACCRETION DISK PARAMETERS (зі статті ApJ Letters 993:L22)
#! ============================================
#! H/r = 0.1 (scale height, як у симуляціях статті)
#! Диск лежить у площині x-z, перпендикулярній до осі y
#! На Figure 1 показано "one scale height" як помаранчевий шахматний візерунок
H_over_r = 0.1  # Scale height H/r = 0.1 (fiducial value from paper)
disk_inner_radius = 3.0  # Inner radius ~3 r_g (близько до горизонту подій)
disk_outer_radius = r0 * 0.8  # Outer radius ~80% від бінарної відстані
disk_color = '#FF8C00'  # Dark orange (як у статті)
disk_alpha = 0.9

#! ============================================
#! PLOT 1: Orbit in x-y plane - "ROSETTE" pattern
#! ============================================
logging.info("Plotting orbital precession rosette (x-y plane)...")

fig, ax = plt.subplots(figsize=(12, 12))

# Color by orbit number to show precession
# Using RELATIVE coordinates (x_rel, y_rel) for proper rosette in all quadrants
norm = plt.Normalize(orbits.min(), orbits.max())
points = np.array([x_rel, y_rel]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

lc = LineCollection(segments, linewidths=0.8, cmap='rainbow', norm=norm)
lc.set_array(orbits)
line = ax.add_collection(lc)
cbar = fig.colorbar(line, ax=ax, label='Orbit number', pad=0.02)

# Set symmetric limits to show all four quadrants
max_extent = max(abs(x_rel.min()), abs(x_rel.max()), abs(y_rel.min()), abs(y_rel.max())) * 1.1
ax.set_xlim(-max_extent, max_extent)
ax.set_ylim(-max_extent, max_extent)
ax.set_xlabel('y/M (relative)', fontsize=14)
ax.set_ylabel('z/M (relative)', fontsize=14)
if orbital_inclination != 0.0:
    ax.set_title(f'OJ 287: Perihelion Precession (Rosette Pattern)\nOrbital Inclination: {orbital_inclination}° around {rotation_axis.upper()}-axis', fontsize=16)
else:
    ax.set_title('OJ 287: Perihelion Precession (Rosette Pattern)\nPure Relativistic 1PN Effect', fontsize=16)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

# Draw accretion disk (edge-on view in x-y plane, disk is in x-z plane)
# Disk extends across the entire x-axis (same length as axis)
disk_thickness_2d = max_extent * 0.01  # Thin disk line

# Left part of disk (from -max_extent to -disk_inner_radius)
ax.fill_between([-max_extent, -disk_inner_radius], [-disk_thickness_2d, -disk_thickness_2d], [disk_thickness_2d, disk_thickness_2d], 
                color=disk_color, alpha=disk_alpha, label='Accretion Disk')
# Right part of disk (from disk_inner_radius to max_extent)
ax.fill_between([disk_inner_radius, max_extent], [-disk_thickness_2d, -disk_thickness_2d], [disk_thickness_2d, disk_thickness_2d], 
                color=disk_color, alpha=disk_alpha)

# Mark the primary black hole at origin
ax.plot(0, 0, 'ko', markersize=15, label='Primary BH (M₁)')
ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.legend(loc='upper right')

fig.tight_layout()
plt.savefig('Precession_Rosette_' + filenameprefix + '.png', dpi=150)
plt.close()
logging.info("Saved: Precession_Rosette_%s.png", filenameprefix)

#! ============================================
#! PLOT 2: Zoomed rosette showing precession clearly
#! ============================================
logging.info("Plotting zoomed precession pattern...")

fig, ax = plt.subplots(figsize=(12, 12))

lc2 = LineCollection(segments, linewidths=1.2, cmap='viridis', norm=norm)
lc2.set_array(orbits)
line2 = ax.add_collection(lc2)
cbar2 = fig.colorbar(line2, ax=ax, label='Orbit number', pad=0.02)

# Zoom to show precession detail (symmetric around origin)
zoom_extent = max_extent * 0.3
ax.set_xlim(-zoom_extent, zoom_extent)
ax.set_ylim(-zoom_extent, zoom_extent)
ax.set_xlabel('y/M (relative)', fontsize=14)
ax.set_ylabel('z/M (relative)', fontsize=14)
if orbital_inclination != 0.0:
    ax.set_title(f'OJ 287: Perihelion Precession (Zoomed)\nInclination: {orbital_inclination}° | {orbitsmax} orbits', fontsize=16)
else:
    ax.set_title('OJ 287: Perihelion Precession (Zoomed)\nShowing orbital shift over %d orbits' % orbitsmax, fontsize=16)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)

# Draw accretion disk (edge-on view) - extends across entire x-axis
zoom_disk_thickness = zoom_extent * 0.01  # Thin disk line

# Left part of disk (from -zoom_extent to -disk_inner_radius)
ax.fill_between([-zoom_extent, -disk_inner_radius], [-zoom_disk_thickness, -zoom_disk_thickness], [zoom_disk_thickness, zoom_disk_thickness], 
                color=disk_color, alpha=disk_alpha, label='Accretion Disk')
# Right part of disk (from disk_inner_radius to zoom_extent)
ax.fill_between([disk_inner_radius, zoom_extent], [-zoom_disk_thickness, -zoom_disk_thickness], [zoom_disk_thickness, zoom_disk_thickness], 
                color=disk_color, alpha=disk_alpha)

ax.plot(0, 0, 'ko', markersize=10, label='Primary BH')
ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.legend(loc='upper right')

fig.tight_layout()
plt.savefig('Precession_Zoomed_' + filenameprefix + '.png', dpi=150)
plt.close()
logging.info("Saved: Precession_Zoomed_%s.png", filenameprefix)

#! ============================================
#! PLOT 3: Perihelion angle vs orbit number
#! ============================================
logging.info("Calculating perihelion angles...")

# Find perihelion positions (minimum r for each orbit)
perihelion_angles = []
perihelion_orbits = []
perihelion_r = []

# Group data by orbit and find minimum r
current_orbit = 0
min_r = float('inf')
min_idx = 0

for i in range(len(r)):
	orbit_num = int(orbits[i])
	if orbit_num > current_orbit:
		if min_r < float('inf'):
			# Calculate perihelion angle using RELATIVE coordinates
			angle = np.arctan2(y_rel[min_idx], x_rel[min_idx]) * 180 / np.pi
			perihelion_angles.append(angle)
			perihelion_orbits.append(current_orbit)
			perihelion_r.append(min_r)
		current_orbit = orbit_num
		min_r = r[i]
		min_idx = i
	elif r[i] < min_r:
		min_r = r[i]
		min_idx = i

# Add last orbit
if min_r < float('inf') and current_orbit > 0:
	angle = np.arctan2(y_rel[min_idx], x_rel[min_idx]) * 180 / np.pi
	perihelion_angles.append(angle)
	perihelion_orbits.append(current_orbit)
	perihelion_r.append(min_r)

# Unwrap angles to show continuous precession
perihelion_angles = np.array(perihelion_angles)
perihelion_angles_unwrapped = np.unwrap(perihelion_angles * np.pi / 180) * 180 / np.pi

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Perihelion angle vs orbit
ax1.plot(perihelion_orbits, perihelion_angles_unwrapped, 'b.-', linewidth=2, markersize=8)
ax1.set_xlabel('Orbit number', fontsize=14)
ax1.set_ylabel('Perihelion angle [degrees]', fontsize=14)
ax1.set_title('OJ 287: Perihelion Angle Evolution', fontsize=16)
ax1.grid(True, alpha=0.3)

# Calculate precession rate
precession_rate = 0
if len(perihelion_orbits) > 1:
	precession_rate = (perihelion_angles_unwrapped[-1] - perihelion_angles_unwrapped[0]) / (perihelion_orbits[-1] - perihelion_orbits[0])
	ax1.text(0.05, 0.95, f'Precession rate: {precession_rate:.2f}°/orbit', 
			transform=ax1.transAxes, fontsize=12, verticalalignment='top',
			bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 2: Perihelion distance vs orbit
ax2.plot(perihelion_orbits, perihelion_r, 'r.-', linewidth=2, markersize=8)
ax2.set_xlabel('Orbit number', fontsize=14)
ax2.set_ylabel('Perihelion distance r/M', fontsize=14)
ax2.set_title('OJ 287: Perihelion Distance (should be constant without radiation)', fontsize=16)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
plt.savefig('Precession_Analysis_' + filenameprefix + '.png', dpi=150)
plt.close()
logging.info("Saved: Precession_Analysis_%s.png", filenameprefix)

#! ============================================
#! PLOT 4: 3D orbit visualization
#! ============================================
logging.info("Plotting 3D orbit...")

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Using RELATIVE coordinates for 3D plot
norm = plt.Normalize(orbits.min(), orbits.max())
points3d = np.array([x_rel, y_rel, z_rel]).T.reshape(-1, 1, 3)
segments3d = np.concatenate([points3d[:-1], points3d[1:]], axis=1)

lc3d = Line3DCollection(segments3d, linewidths=0.8, cmap='rainbow')
lc3d.set_array(orbits)
line3d = ax.add_collection(lc3d)
fig.colorbar(line3d, ax=ax, label='Orbit number', pad=0.1)

# Draw accretion disk in 3D (disk in x-z plane, perpendicular to y-axis)
# Create annulus (ring) in x-z plane
theta_disk = np.linspace(0, 2*np.pi, 100)
# Outer edge of disk
x_disk_outer = disk_outer_radius * np.cos(theta_disk)
z_disk_outer = disk_outer_radius * np.sin(theta_disk)
y_disk_outer = np.zeros_like(theta_disk)
# Inner edge of disk
x_disk_inner = disk_inner_radius * np.cos(theta_disk)
z_disk_inner = disk_inner_radius * np.sin(theta_disk)
y_disk_inner = np.zeros_like(theta_disk)

# Plot disk as surface (annulus in x-z plane)
from matplotlib.colors import LinearSegmentedColormap

# Create radial points for disk surface
n_radial = 20
radii = np.linspace(disk_inner_radius, disk_outer_radius, n_radial)
theta_surf = np.linspace(0, 2*np.pi, 50)
R, Theta = np.meshgrid(radii, theta_surf)
X_disk = R * np.cos(Theta)
Z_disk = R * np.sin(Theta)
Y_disk = np.zeros_like(X_disk)

# Custom colormap for accretion disk (hot colors)
disk_cmap = LinearSegmentedColormap.from_list('accretion', 
    ['darkred', 'red', 'orange', 'yellow', 'white'], N=256)

# Color by radius (hotter near center)
disk_colors = (R.max() - R) / (R.max() - R.min())
ax.plot_surface(X_disk, Y_disk, Z_disk, facecolors=disk_cmap(disk_colors), 
                alpha=0.7, linewidth=0, antialiased=True, shade=False)

# Add disk outline
ax.plot(x_disk_outer, y_disk_outer, z_disk_outer, color='orange', linewidth=1.5, alpha=0.8)
ax.plot(x_disk_inner, y_disk_inner, z_disk_inner, color='yellow', linewidth=1.5, alpha=0.8)

# Symmetric limits for all quadrants
max_xy = max(abs(x_rel.min()), abs(x_rel.max()), abs(y_rel.min()), abs(y_rel.max()))
max_z = max(abs(z_rel.min()), abs(z_rel.max())) if abs(z_rel).max() > 0 else max_xy * 0.1
ax.set_xlim(-max_xy, max_xy)
ax.set_ylim(-max_xy, max_xy)
ax.set_zlim(-max_z, max_z)

ax.set_xlabel('x/M (relative)')
ax.set_ylabel('y/M (relative)')
ax.set_zlabel('z/M (relative)')
if orbital_inclination != 0.0:
    ax.set_title(f'OJ 287: 3D Orbital Precession\nInclination: {orbital_inclination}° around {rotation_axis.upper()}-axis')
else:
    ax.set_title('OJ 287: 3D Orbital Precession with Accretion Disk')

# Mark primary at origin
ax.scatter([0], [0], [0], color='black', s=100, label='Primary BH')
ax.legend()

fig.tight_layout()
plt.savefig('Precession_3D_' + filenameprefix + '.png', dpi=150)
plt.close()
logging.info("Saved: Precession_3D_%s.png", filenameprefix)

#! ============================================
#! Summary output
#! ============================================
logging.info("=" * 60)
logging.info("SIMULATION SUMMARY")
logging.info("=" * 60)
logging.info("Total orbits simulated: %.1f", orbits[-1])
logging.info("Total simulation time: %.2f years", ty[-1])
if orbital_inclination != 0.0:
    logging.info("Orbital inclination: %.1f degrees around %s-axis", orbital_inclination, rotation_axis.upper())
if len(perihelion_orbits) > 1:
	total_precession = perihelion_angles_unwrapped[-1] - perihelion_angles_unwrapped[0]
	logging.info("Total perihelion precession: %.2f degrees", total_precession)
	logging.info("Precession rate: %.2f degrees/orbit", precession_rate)
	logging.info("Precession rate: %.2f degrees/year", precession_rate * (orbits[-1] / ty[-1]))
logging.info("=" * 60)

#! Move files to data directory
now = datetime.datetime.now()
dirname = 'data/' + filenameprefix + '_' + now.strftime("%y-%m-%d-%H%M")

if not os.path.isdir('data'):
	call(['mkdir', 'data'])

call(['mkdir', dirname])
logging.info(dirname + " has been created!")

if os.path.isdir(dirname):
	os.system('mv *' + filenameprefix + '*.png ' + dirname)
	os.system('mv *' + filenameprefix + '*.ini ' + dirname)
	os.system('mv *' + filenameprefix + '*.txt ' + dirname)  # координати
	
	if removedat:
		os.system('rm *' + filenameprefix + '*.dat')
		logging.info("Generated data files has been removed!")
	else:
		os.system('mv *' + filenameprefix + '*.dat ' + dirname)
		logging.info('Generated files moved to ' + dirname)
else:
	logging.warning("Directory doesn't exist!")

logging.info("Done! Check the 'data' folder for results.")
