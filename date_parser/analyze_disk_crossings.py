#!/usr/bin/env python3
"""
OJ 287 Binary Black Hole - Accretion Disk Crossing Events Analyzer

This script parses coordinate data and identifies when the secondary black hole
crosses the accretion disk plane (z_rel = 0) of the primary black hole.

A disk crossing occurs when z_rel changes sign.
- Upward crossing (z: - → +): z_rel goes from negative to positive
- Downward crossing (z: + → -): z_rel goes from positive to negative

Time Scaling:
  The simulation may use arbitrary time units. For OJ287, we rescale time
  based on the known orbital period of ~12 years to get physically correct results.
"""

import os
import sys
from dataclasses import dataclass
from typing import List, Tuple

# OJ287 physical parameters
OJ287_ORBITAL_PERIOD_YEARS = 12.0  # Known orbital period of OJ287 in years


@dataclass
class DiskCrossing:
    """Represents a single disk crossing event."""
    number: int
    time_before: float
    time_after: float
    direction: str  # "up" or "down"
    separation: float
    orbit_number: float


def parse_coordinates_file(filepath: str, use_physical_time: bool = True) -> Tuple[List[Tuple[float, float, float, float]], dict]:
    """
    Parse the coordinates file and extract relevant data.
    
    Args:
        filepath: Path to the coordinates file
        use_physical_time: If True, rescale time based on OJ287's 12-year orbital period
    
    Returns:
        - List of tuples: (time, z_rel, separation, orbit_number)
        - Dictionary with metadata
    """
    data = []
    metadata = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Parse header comments for metadata
            if line.startswith('#'):
                if 'M1 =' in line:
                    metadata['M1'] = line.split('=')[1].split()[0].strip()
                elif 'M2 =' in line:
                    metadata['M2'] = line.split('=')[1].split()[0].strip()
                continue
            
            # Parse data line
            parts = line.split()
            if len(parts) >= 12:
                try:
                    sim_time = float(parts[0])
                    z_rel = float(parts[8])  # Column 9: z_rel/M
                    separation = float(parts[10])  # Column 11: separation/M
                    orbit_number = float(parts[11])  # Column 12: orbit_number
                    
                    # Calculate physical time based on orbit number and known period
                    if use_physical_time:
                        time = orbit_number * OJ287_ORBITAL_PERIOD_YEARS
                    else:
                        time = sim_time
                    
                    data.append((time, z_rel, separation, orbit_number))
                except (ValueError, IndexError):
                    continue
    
    return data, metadata


def find_disk_crossings(data: List[Tuple[float, float, float, float]]) -> List[DiskCrossing]:
    """
    Find all disk crossing events where z_rel changes sign.
    
    Args:
        data: List of (time, z_rel, separation, orbit_number) tuples
    
    Returns:
        List of DiskCrossing objects
    """
    crossings = []
    crossing_number = 0
    
    for i in range(1, len(data)):
        prev_time, prev_z, prev_sep, prev_orbit = data[i-1]
        curr_time, curr_z, curr_sep, curr_orbit = data[i]
        
        # Check for sign change (crossing)
        if prev_z * curr_z < 0:  # Signs are different
            crossing_number += 1
            
            # Determine direction
            if prev_z < 0 and curr_z > 0:
                direction = "up"
            else:
                direction = "down"
            
            # Use values closer to the crossing point (interpolate separation and orbit)
            avg_separation = (prev_sep + curr_sep) / 2
            avg_orbit = (prev_orbit + curr_orbit) / 2
            
            crossing = DiskCrossing(
                number=crossing_number,
                time_before=prev_time,
                time_after=curr_time,
                direction=direction,
                separation=avg_separation,
                orbit_number=avg_orbit
            )
            crossings.append(crossing)
    
    return crossings


def calculate_intervals(crossings: List[DiskCrossing]) -> List[Tuple[int, int, float]]:
    """
    Calculate time intervals between consecutive crossings.
    
    Returns:
        List of (crossing1, crossing2, interval_in_years) tuples
    """
    intervals = []
    
    for i in range(1, len(crossings)):
        prev_crossing = crossings[i-1]
        curr_crossing = crossings[i]
        
        # Use the midpoint time of each crossing
        prev_time = (prev_crossing.time_before + prev_crossing.time_after) / 2
        curr_time = (curr_crossing.time_before + curr_crossing.time_after) / 2
        
        interval = curr_time - prev_time
        intervals.append((prev_crossing.number, curr_crossing.number, interval))
    
    return intervals


def generate_report(filepath: str, output_path: str = None, use_physical_time: bool = True):
    """
    Generate the disk crossing analysis report.
    
    Args:
        filepath: Path to the coordinates file
        output_path: Optional path for output file. If None, prints to stdout.
        use_physical_time: If True, use OJ287's 12-year orbital period for time scaling
    """
    # Parse data
    data, metadata = parse_coordinates_file(filepath, use_physical_time)
    
    if not data:
        print("Error: No data found in file.")
        return
    
    # Find crossings
    crossings = find_disk_crossings(data)
    
    # Calculate statistics
    total_points = len(data)
    time_min = data[0][0]
    time_max = data[-1][0]
    total_crossings = len(crossings)
    upward_crossings = sum(1 for c in crossings if c.direction == "up")
    downward_crossings = sum(1 for c in crossings if c.direction == "down")
    
    # Calculate intervals
    intervals = calculate_intervals(crossings)
    
    # Get orbit range
    orbit_min = data[0][3]
    orbit_max = data[-1][3]
    
    # Build report
    lines = []
    lines.append("")
    lines.append("OJ 287 Binary Black Hole - Accretion Disk Crossing Events")
    lines.append("=" * 60)
    lines.append("")
    lines.append(f"Data points analyzed: {total_points}")
    lines.append(f"Time span: {time_min:.4f} - {time_max:.4f} years")
    lines.append(f"Orbital period: {OJ287_ORBITAL_PERIOD_YEARS:.1f} years (OJ287 known value)")
    lines.append(f"Orbits covered: {orbit_min:.4f} - {orbit_max:.4f}")
    lines.append(f"Total disk crossings found: {total_crossings}")
    lines.append(f"  - Upward crossings (z: - → +): {upward_crossings}")
    lines.append(f"  - Downward crossings (z: + → -): {downward_crossings}")
    lines.append("")
    lines.append("-" * 60)
    lines.append("Disk Crossing Events (chronological order):")
    lines.append("-" * 60)
    lines.append("")
    
    # Table header
    lines.append(f"{'#':>4}  {'Time Before':>12}  {'Time After':>12}  {'Direction':>10}  {'Separation':>12}  {'Orbit #':>10}")
    lines.append("")
    
    # Table rows
    for c in crossings:
        lines.append(
            f"{c.number:>4}  {c.time_before:>12.6f}  {c.time_after:>12.6f}  "
            f"{c.direction:>10}  {c.separation:>12.4f}  {c.orbit_number:>10.4f}"
        )
    
    lines.append("")
    lines.append("-" * 60)
    lines.append("Time intervals between consecutive crossings:")
    lines.append("-" * 60)
    lines.append("")
    
    for prev_num, curr_num, interval in intervals:
        lines.append(f"Crossing {prev_num:>2} -> {curr_num:>2}: {interval:.6f} years")
    
    lines.append("")
    
    # Output report
    report = "\n".join(lines)
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(report)
        print(f"Report saved to: {output_path}")
    else:
        print(report)


def find_latest_data_folder(base_dir: str) -> str:
    """
    Find the latest OJ287_precession_* folder in the data directory.
    
    Args:
        base_dir: Path to the data directory
    
    Returns:
        Path to the latest folder, or None if not found
    """
    import glob
    
    # Find all matching folders
    pattern = os.path.join(base_dir, "OJ287_precession_*")
    folders = glob.glob(pattern)
    
    if not folders:
        return None
    
    # Sort by modification time (newest first)
    folders.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    
    return folders[0]


def main():
    """Main entry point."""
    # Base data directory
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(project_root, "data")
    
    # Parse command line arguments
    if len(sys.argv) >= 2:
        input_file = sys.argv[1]
    else:
        # Find the latest data folder automatically
        latest_folder = find_latest_data_folder(data_dir)
        
        if latest_folder is None:
            print(f"Error: No OJ287_precession_* folders found in {data_dir}")
            sys.exit(1)
        
        input_file = os.path.join(latest_folder, "OJ287_precession_coordinates.txt")
        print(f"Using latest data folder: {os.path.basename(latest_folder)}")
    
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
    else:
        # Default output to same directory as input
        input_dir = os.path.dirname(input_file)
        output_file = os.path.join(input_dir, "results.txt")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
    
    # Generate report
    generate_report(input_file, output_file)


if __name__ == "__main__":
    main()
