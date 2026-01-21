#!/usr/bin/env python3
"""
Parse Disk Crossings from CBwaves DAT file

This script reads the output from CBwaves simulation (.dat file)
and finds all moments when the secondary black hole crosses through
the accretion disk of the primary black hole.

Scientific approach for OJ 287-like systems:
- Orbital period ~12 years
- Two disk crossings per orbit (double flare)
- Crossings are grouped by orbital period to filter noise

Usage:
    python parse_disk_crossings.py input.dat [output.txt]

If output file is not specified, it will be auto-generated based on input filename.
"""

import os
import sys
import csv
import logging
import argparse

# Logging setup
logging.basicConfig(
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO
)

# Physical constants
YEAR_SECONDS = 365.25636 * 24.0 * 3600.0  # Year in seconds

# OJ 287 orbital parameters
ORBITAL_PERIOD = 12.0  # Orbital period in years
HALF_PERIOD = 6.0  # Half orbital period (time between disk passages)
DOUBLE_CROSSING_WINDOW = 2.0  # Max time between entry/exit crossings (years)
MIN_PASSAGE_GAP = 3.0  # Minimum time between disk passages (years)


def parse_disk_crossings(dat_file, txt_file):
    """
    Read coordinates from DAT file and find disk crossings.
    
    DAT file format (columns):
    - col 0: time
    - col 1: x1 (primary BH x-coordinate)
    - col 2: y1 (primary BH y-coordinate)
    - col 3: z1 (primary BH z-coordinate)
    - col 4: x2 (secondary BH x-coordinate)
    - col 5: y2 (secondary BH y-coordinate) <- used for disk crossing detection
    - col 6: z2 (secondary BH z-coordinate)
    
    Parameters:
    -----------
    dat_file : str
        Path to the input .dat file from CBwaves simulation
    txt_file : str
        Path to the output .txt file for crossings
        
    Returns:
    --------
    int : Number of crossings found, or -1 on error
    """
    
    logging.info("=" * 60)
    logging.info("PARSING DISK CROSSINGS")
    logging.info("=" * 60)
    logging.info("Input file: %s", dat_file)
    logging.info("Output file: %s", txt_file)
    
    # Check if input file exists
    if not os.path.isfile(dat_file):
        logging.error("Input file not found: %s", dat_file)
        return -1
    
    # ============================================
    # Read data from DAT file
    # ============================================
    logging.info("Reading data from DAT file...")
    
    tt = []   # time
    yy1 = []  # y-coordinate of primary BH
    yy2 = []  # y-coordinate of secondary BH
    
    with open(dat_file) as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            try:
                # Filter out empty strings from row
                row = [x for x in row if x]
                tt.append(float(row[0]))
                yy1.append(float(row[2]))  # y1 coordinate (col 2)
                yy2.append(float(row[5]))  # y2 coordinate (col 5)
            except (ValueError, IndexError):
                # Skip malformed lines
                continue
    
    total_points = len(tt)
    logging.info("Read %d data points from DAT file", total_points)
    
    if total_points == 0:
        logging.error("No data points read from DAT file!")
        return -1
    
    # Convert time to model years
    ty = [t / YEAR_SECONDS for t in tt]
    
    # Calculate RELATIVE y-coordinate (secondary relative to primary)
    # This is scientifically correct for disk crossing detection
    y_rel = [yy2[i] - yy1[i] for i in range(len(yy2))]
    
    # ============================================
    # Find ALL zero-crossings of y_rel
    # ============================================
    logging.info("Searching for disk crossings (y_rel crossing 0)...")
    
    raw_crossings = []
    
    for i in range(len(y_rel) - 1):
        # Check if y_rel changes sign (crosses 0)
        if y_rel[i] * y_rel[i + 1] < 0:
            # Determine direction
            if y_rel[i] > 0 and y_rel[i + 1] < 0:
                direction = "down"
            else:
                direction = "up"
            
            raw_crossings.append({
                'time_before': ty[i],
                'time_after': ty[i + 1],
                'y_rel_before': y_rel[i],
                'y_rel_after': y_rel[i + 1],
                'direction': direction
            })
    
    logging.info("Raw crossings found: %d", len(raw_crossings))
    
    if len(raw_crossings) == 0:
        logging.warning("No disk crossings found in the data!")
        return 0
    
    # ============================================
    # STEP 1: Group crossings into disk passages
    # Each passage has 1-2 crossings (entry/exit) within ~2 years
    # ============================================
    logging.info("Applying orbital period filter (P = %.1f years)...", ORBITAL_PERIOD)
    
    passages = []
    current_passage = [raw_crossings[0]]
    
    for i in range(1, len(raw_crossings)):
        time_gap = raw_crossings[i]['time_before'] - current_passage[-1]['time_before']
        
        if time_gap < DOUBLE_CROSSING_WINDOW:
            # Same disk passage - add to current group
            current_passage.append(raw_crossings[i])
        else:
            # New disk passage - save current and start new
            passages.append(current_passage)
            current_passage = [raw_crossings[i]]
    
    passages.append(current_passage)
    
    logging.info("Disk passages identified: %d", len(passages))
    
    # ============================================
    # STEP 2: Group passages into full orbits (~12 years)
    # Each orbit has 2 passages separated by ~6 years
    # ============================================
    orbits = []
    current_orbit = [passages[0]]
    
    for i in range(1, len(passages)):
        # Time from start of current orbit to this passage
        time_from_orbit_start = passages[i][0]['time_before'] - current_orbit[0][0]['time_before']
        
        if time_from_orbit_start < ORBITAL_PERIOD - 1:  # Within same orbit (with 1-year tolerance)
            current_orbit.append(passages[i])
        else:
            # New orbit
            orbits.append(current_orbit)
            current_orbit = [passages[i]]
    
    orbits.append(current_orbit)
    
    logging.info("Full orbits identified: %d", len(orbits))
    
    # ============================================
    # STEP 3: Extract all crossings with orbit/passage info
    # ============================================
    crossings = []
    
    for orbit_num, orbit_passages in enumerate(orbits, 1):
        for passage_num, passage in enumerate(orbit_passages, 1):
            passage_label = "A" if passage_num == 1 else "B"
            
            if len(passage) >= 2:
                # Double crossing (entry + exit)
                entry = passage[0]
                exit_cr = passage[-1]
                
                crossing_id = len(crossings) + 1
                crossings.append({
                    'id': crossing_id,
                    'orbit': orbit_num,
                    'passage': passage_label,
                    'type': 'entry',
                    'time_before': entry['time_before'],
                    'time_after': entry['time_after'],
                    'y_rel_before': entry['y_rel_before'],
                    'y_rel_after': entry['y_rel_after'],
                    'direction': entry['direction']
                })
                
                crossing_id = len(crossings) + 1
                crossings.append({
                    'id': crossing_id,
                    'orbit': orbit_num,
                    'passage': passage_label,
                    'type': 'exit',
                    'time_before': exit_cr['time_before'],
                    'time_after': exit_cr['time_after'],
                    'y_rel_before': exit_cr['y_rel_before'],
                    'y_rel_after': exit_cr['y_rel_after'],
                    'direction': exit_cr['direction']
                })
                
                interval = exit_cr['time_before'] - entry['time_before']
                logging.info(
                    "  Orbit %d, Passage %s: entry=%.2f, exit=%.2f yrs (interval=%.2f yrs)",
                    orbit_num, passage_label, entry['time_before'], exit_cr['time_before'], interval
                )
            else:
                # Single crossing
                single = passage[0]
                crossing_id = len(crossings) + 1
                crossings.append({
                    'id': crossing_id,
                    'orbit': orbit_num,
                    'passage': passage_label,
                    'type': 'single',
                    'time_before': single['time_before'],
                    'time_after': single['time_after'],
                    'y_rel_before': single['y_rel_before'],
                    'y_rel_after': single['y_rel_after'],
                    'direction': single['direction']
                })
                
                logging.info(
                    "  Orbit %d, Passage %s: single crossing at %.2f yrs",
                    orbit_num, passage_label, single['time_before']
                )
    
    num_crossings = len(crossings)
    num_orbits = len(orbits)
    num_passages = len(passages)
    
    logging.info("-" * 60)
    logging.info("SCIENTIFIC SUMMARY:")
    logging.info("  Full orbits (P=%.1f yrs): %d", ORBITAL_PERIOD, num_orbits)
    logging.info("  Disk passages (2 per orbit): %d", num_passages)
    logging.info("  Total crossings (2 per passage): %d", num_crossings)
    
    # Count directions
    up_count = sum(1 for c in crossings if c['direction'] == 'up')
    down_count = sum(1 for c in crossings if c['direction'] == 'down')
    logging.info("  Upward crossings: %d", up_count)
    logging.info("  Downward crossings: %d", down_count)
    
    # Verify orbital period
    if num_orbits >= 2:
        orbit_start_times = [orbit[0][0]['time_before'] for orbit in orbits]
        periods = [orbit_start_times[i+1] - orbit_start_times[i] for i in range(len(orbit_start_times)-1)]
        avg_period = sum(periods) / len(periods)
        logging.info("  Measured orbital period: %.2f years", avg_period)
        logging.info("  Expected period: %.2f years", ORBITAL_PERIOD)
        
        # Verify half-period (between passages within same orbit)
        half_periods = []
        for orbit in orbits:
            if len(orbit) >= 2:
                hp = orbit[1][0]['time_before'] - orbit[0][0]['time_before']
                half_periods.append(hp)
        if half_periods:
            avg_half = sum(half_periods) / len(half_periods)
            logging.info("  Measured half-period: %.2f years", avg_half)
            logging.info("  Expected half-period: %.2f years", HALF_PERIOD)
    
    # ============================================
    # Write crossings to output file
    # ============================================
    logging.info("Writing crossings to output file...")
    
    with open(txt_file, 'w') as f:
        # Header
        f.write("# Disk Crossings - Secondary BH through Primary BH Accretion Disk\n")
        f.write("# Generated by parse_disk_crossings.py\n")
        f.write("#\n")
        f.write("# Scientific Parameters (OJ 287 system):\n")
        f.write("#   Orbital period: %.1f years\n" % ORBITAL_PERIOD)
        f.write("#   Half-period (between passages): %.1f years\n" % HALF_PERIOD)
        f.write("#   Double crossing window: %.1f years\n" % DOUBLE_CROSSING_WINDOW)
        f.write("#\n")
        f.write("# Input file: %s\n" % os.path.basename(dat_file))
        f.write("# Model duration: %.2f years\n" % ty[-1])
        f.write("# Total data points: %d\n" % total_points)
        f.write("#\n")
        f.write("# Results:\n")
        f.write("#   Raw crossings detected: %d\n" % len(raw_crossings))
        f.write("#   Full orbits: %d\n" % num_orbits)
        f.write("#   Disk passages: %d (2 per orbit)\n" % num_passages)
        f.write("#   Significant crossings: %d (2 per passage)\n" % num_crossings)
        f.write("#   Upward crossings: %d\n" % up_count)
        f.write("#   Downward crossings: %d\n" % down_count)
        f.write("#\n")
        f.write("# Method: Relative y-coordinate (y_rel = y2 - y1) crossing zero\n")
        f.write("# Structure: Orbit -> Passage (A or B) -> Crossing (entry/exit)\n")
        f.write("#\n")
        f.write("# Column format (9 columns):\n")
        f.write("# 1: crossing_id - Sequential crossing number\n")
        f.write("# 2: orbit - Full orbital cycle number\n")
        f.write("# 3: passage - Disk passage within orbit (A=first, B=second)\n")
        f.write("# 4: type - 'entry', 'exit', or 'single'\n")
        f.write("# 5: time_before[years] - Model time just before crossing\n")
        f.write("# 6: time_after[years] - Model time just after crossing\n")
        f.write("# 7: direction - 'up' (y_rel: - -> +) or 'down' (y_rel: + -> -)\n")
        f.write("# 8: y_rel_before - Relative y-coordinate before crossing\n")
        f.write("# 9: y_rel_after - Relative y-coordinate after crossing\n")
        f.write("#\n")
        f.write("#" + "=" * 120 + "\n")
        
        # Data
        for c in crossings:
            f.write("%d  %d  %s  %s  %.4f  %.4f  %s  %.10e  %.10e\n" % (
                c['id'],
                c['orbit'],
                c['passage'],
                c['type'],
                c['time_before'],
                c['time_after'],
                c['direction'],
                c['y_rel_before'],
                c['y_rel_after']
            ))
    
    logging.info("Crossings saved to: %s", txt_file)
    logging.info("Output file size: %.2f KB", os.path.getsize(txt_file) / 1024)
    
    # ============================================
    # Summary
    # ============================================
    logging.info("=" * 60)
    logging.info("FINAL SUMMARY")
    logging.info("=" * 60)
    logging.info("Model time range: %.2f to %.2f years", ty[0], ty[-1])
    logging.info("Full orbits: %d", num_orbits)
    logging.info("Disk passages: %d", num_passages)
    logging.info("Total crossings: %d", num_crossings)
    if crossings:
        logging.info("First crossing: %.2f years (Orbit %d, Passage %s)", 
                     crossings[0]['time_before'], crossings[0]['orbit'], crossings[0]['passage'])
        logging.info("Last crossing: %.2f years (Orbit %d, Passage %s)", 
                     crossings[-1]['time_before'], crossings[-1]['orbit'], crossings[-1]['passage'])
    logging.info("=" * 60)
    
    return num_crossings


def main():
    """Main entry point with command line argument parsing."""
    
    parser = argparse.ArgumentParser(
        description='Parse disk crossings from CBwaves DAT file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python parse_disk_crossings.py OJ287_precession.dat
    python parse_disk_crossings.py OJ287_precession.dat crossings.txt
    python parse_disk_crossings.py data/folder/simulation.dat
        """
    )
    
    parser.add_argument('input_dat',
                        help='Input DAT file from CBwaves simulation')
    parser.add_argument('output_txt', nargs='?', default=None,
                        help='Output TXT file for crossings (auto-generated if not specified)')
    
    args = parser.parse_args()
    
    # Auto-generate output filename if not specified
    if args.output_txt is None:
        base_name = os.path.splitext(args.input_dat)[0]
        args.output_txt = base_name + '_crossings.txt'
    
    result = parse_disk_crossings(args.input_dat, args.output_txt)
    
    if result >= 0:
        logging.info("Done! Found %d disk crossings.", result)
    else:
        logging.error("Parsing failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
