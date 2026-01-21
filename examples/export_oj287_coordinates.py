#!/usr/bin/env python3
"""
Export OJ 287 Coordinates from DAT to TXT

This script reads the output from CBwaves simulation (.dat file)
and exports the coordinates to a text file (.txt) for further processing
and visualization by oj287_procession.py.

Usage:
    python export_oj287_coordinates.py [input.dat] [output.txt]

If no arguments provided, uses default filenames.
"""

#* Mathematical imports
import math
import numpy as np

#* Utility imports
import os
import sys
import csv
import logging
import argparse

#* Logging setup
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

#* Definition of auxiliary variables
pi = 3.141592653589793      #* Pi
c = 2.99792458e8            #* Speed of light in [m/s]
G = 6.674184e-11            #* Gravitational constant 
Gc2 = G/(c*c)               #* G/c^2
ly = 9460730472580800.0     #* Light year
pc = 3.26156 * ly           #* Parsec
year = 365.25636*24.*3600.  #* Orbital period time of Earth in seconds
msun = 1476.62504           #* Mass of the sun in G=c=1 system, [m]


def export_coordinates(dat_file, txt_file, max_points=10000):
    """
    Read coordinates from DAT file and export to TXT file.
    
    Parameters:
    -----------
    dat_file : str
        Path to the input .dat file from CBwaves simulation
    txt_file : str
        Path to the output .txt file for coordinates
    max_points : int
        Maximum number of data points to export (for file size control)
        
    Returns:
    --------
    str : Path to the created txt file
    """
    
    #! ============================================
    #! ПАРАМЕТРИ ЗІ СТАТТІ (must match oj287_procession.py)
    #! q = 0.1, r0 = 16.4 r_g
    #! ============================================
    M1 = 1.8348*math.pow(10, 10)  # Primary SMBH in solar masses
    M2 = M1 * 0.1                 # Secondary BH: q = M2/M1 = 0.1
    m1 = M1*msun                   # in meters (G=c=1)
    m2 = M2*msun
    M = m1 + m2                    # Total mass
    
    logging.info("=" * 60)
    logging.info("EXPORTING OJ 287 COORDINATES")
    logging.info("=" * 60)
    logging.info("Input file: %s", dat_file)
    logging.info("Output file: %s", txt_file)
    
    # Check if input file exists
    if not os.path.isfile(dat_file):
        logging.error("Input file not found: %s", dat_file)
        return None
    
    #! ============================================
    #! Read data from DAT file
    #! ============================================
    logging.info("Reading data from DAT file...")
    
    tt = []
    xx1 = []  # coordinates of primary BH
    yy1 = []
    zz1 = []
    xx2 = []  # coordinates of secondary BH
    yy2 = []
    zz2 = []
    rr = []
    orbits_list = []
    ecc_list = []
    
    with open(dat_file) as zw:
        zwval = csv.reader(zw, delimiter=' ')
        for row in zwval:
            try:
                tt.append(float(row[0]))
                xx1.append(float(row[1]))  # x1 coordinate of primary BH
                yy1.append(float(row[2]))  # y1 coordinate of primary BH
                zz1.append(float(row[3]))  # z1 coordinate of primary BH
                xx2.append(float(row[4]))  # x2 coordinate of secondary BH
                yy2.append(float(row[5]))  # y2 coordinate of secondary BH
                zz2.append(float(row[6]))  # z2 coordinate of secondary BH
                rr.append(float(row[7]))   # separation r
                # row[8] is ecc_r
                if len(row) > 8:
                    ecc_list.append(float(row[8]))
                else:
                    ecc_list.append(0.0)
                orbits_list.append(float(row[9]))  # orbit number
            except (ValueError, IndexError) as e:
                # Skip malformed lines
                continue
    
    logging.info("Read %d data points from DAT file", len(tt))
    
    if len(tt) == 0:
        logging.error("No data points read from DAT file!")
        return None
    
    #! ============================================
    #! Normalize coordinates
    #! ============================================
    x1 = np.array([i/M for i in xx1])
    y1 = np.array([j/M for j in yy1])
    z1 = np.array([k/M for k in zz1])
    x2 = np.array([i/M for i in xx2])
    y2 = np.array([j/M for j in yy2])
    z2 = np.array([k/M for k in zz2])
    r = np.array([i/M for i in rr])
    ty = np.array([i/year for i in tt])
    orbits = np.array(orbits_list)
    ecc = np.array(ecc_list)
    
    # Calculate RELATIVE coordinates (secondary relative to primary)
    x_rel = x2 - x1
    y_rel = y2 - y1
    z_rel = z2 - z1
    
    #! ============================================
    #! Write coordinates to TXT file
    #! ============================================
    logging.info("Writing coordinates to TXT file...")
    
    # Calculate step for thinning data if necessary
    step = max(1, len(ty) // max_points)
    num_points = len(range(0, len(ty), step))
    
    with open(txt_file, 'w') as f:
        # Header with metadata
        f.write("# OJ 287 Binary Black Hole Coordinates\n")
        f.write("# Generated by export_oj287_coordinates.py\n")
        f.write("#\n")
        f.write("# System Parameters:\n")
        f.write("# M1 = %.6e solar masses (Primary SMBH)\n" % M1)
        f.write("# M2 = %.6e solar masses (Secondary BH)\n" % M2)
        f.write("# M_total = %.6e meters (G=c=1 units)\n" % M)
        f.write("#\n")
        f.write("# Coordinate normalization: All spatial coordinates divided by M_total\n")
        f.write("# Time units: years\n")
        f.write("#\n")
        f.write("# x_rel, y_rel, z_rel = relative coordinates (secondary - primary)\n")
        f.write("#   These are used for rosette/precession plots centered on primary BH\n")
        f.write("#\n")
        f.write("# Column format (12 columns):\n")
        f.write("# 1: time[years]\n")
        f.write("# 2: y1/M (primary BH y-coordinate)\n")
        f.write("# 3: z1/M (primary BH z-coordinate)\n")
        f.write("# 4: x1/M (primary BH x-coordinate)\n")
        f.write("# 5: y2/M (secondary BH y-coordinate)\n")
        f.write("# 6: z2/M (secondary BH z-coordinate)\n")
        f.write("# 7: x2/M (secondary BH x-coordinate)\n")
        f.write("# 8: y_rel/M (relative y = y2 - y1)\n")
        f.write("# 9: z_rel/M (relative z = z2 - z1)\n")
        f.write("# 10: x_rel/M (relative x = x2 - x1)\n")
        f.write("# 11: separation/M (distance between BHs)\n")
        f.write("# 12: orbit_number\n")
        f.write("#\n")
        f.write("# Data points: %d (out of %d total, step=%d)\n" % (num_points, len(ty), step))
        f.write("#" + "="*120 + "\n")
        
        # Data (with thinning for file size control)
        for i in range(0, len(ty), step):
            f.write("%.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.10e  %.4f\n" % 
                    (ty[i], x1[i], y1[i], z1[i], x2[i], y2[i], z2[i], 
                     x_rel[i], y_rel[i], z_rel[i], r[i], orbits[i]))
    
    logging.info("Coordinates saved: %d data points", num_points)
    logging.info("Output file size: %.2f KB", os.path.getsize(txt_file) / 1024)
    
    #! ============================================
    #! Summary statistics
    #! ============================================
    logging.info("=" * 60)
    logging.info("DATA SUMMARY")
    logging.info("=" * 60)
    logging.info("Time range: %.4f to %.4f years", ty[0], ty[-1])
    logging.info("Orbits range: %.2f to %.2f", orbits[0], orbits[-1])
    logging.info("Separation range: %.4f to %.4f M", r.min(), r.max())
    logging.info("x_rel range: %.4f to %.4f M", x_rel.min(), x_rel.max())
    logging.info("y_rel range: %.4f to %.4f M", y_rel.min(), y_rel.max())
    logging.info("z_rel range: %.4f to %.4f M", z_rel.min(), z_rel.max())
    logging.info("=" * 60)
    
    return txt_file


def main():
    """Main entry point with command line argument parsing."""
    
    parser = argparse.ArgumentParser(
        description='Export OJ 287 coordinates from DAT to TXT format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python export_oj287_coordinates.py
    python export_oj287_coordinates.py OJ287_precession.dat coordinates.txt
    python export_oj287_coordinates.py --max-points 5000 input.dat output.txt
        """
    )
    
    parser.add_argument('input_dat', nargs='?', default='OJ287_precession.dat',
                        help='Input DAT file from CBwaves (default: OJ287_precession.dat)')
    parser.add_argument('output_txt', nargs='?', default='OJ287_precession_coordinates.txt',
                        help='Output TXT file for coordinates (default: OJ287_precession_coordinates.txt)')
    parser.add_argument('--max-points', type=int, default=10000,
                        help='Maximum number of data points to export (default: 10000)')
    
    args = parser.parse_args()
    
    result = export_coordinates(args.input_dat, args.output_txt, args.max_points)
    
    if result:
        logging.info("Done! Coordinates exported to: %s", args.output_txt)
    else:
        logging.error("Export failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
