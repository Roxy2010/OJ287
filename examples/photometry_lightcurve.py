#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OJ 287 Light Curve Plotter with Precession Model Overlay
=========================================================
Plots observed photometric data and overlays predicted outburst times
from the precession model (disk crossing events from results.txt).

Usage:
    python plot_oj287_lightcurve.py

The script will prompt for:
    1. Path to TXT file with photometric data
    2. Filter name (I, V, R, etc.) for the plot title

File format:
    Column 1: Julian Date (JD) - horizontal axis
    Column 2: Magnitude (orange color)
    Column 3: Magnitude (light blue color)
"""

#* Mathematical imports
import math
import numpy as np
import re

#* Utility imports
import os
import sys
import logging
from pathlib import Path

#* Matplotlib imports
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

#* Configure logging
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

#* Matplotlib settings
plt.rc('axes', labelsize=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

#! ============================================
#! REFERENCE EPOCH - First observed outburst peak (t₀)
#! ============================================
T0_JD = 2454145.455  # Reference epoch: first disk crossing (2007.0)

#! ============================================
#! PATH TO RESULTS FILE FROM PRECESSION MODEL
#! ============================================
# This file contains disk crossing times from oj287_procession.py simulation
RESULTS_FILE = None  # Will be auto-detected or specified


def jd_to_calendar(jd):
    """Convert JD to approximate calendar year."""
    return 2000.0 + (jd - 2451545.0) / 365.25


def calendar_to_jd(year_val):
    """Convert calendar year to JD."""
    return 2451545.0 + (year_val - 2000.0) * 365.25


def find_results_file():
    """Find the most recent results.txt file from precession simulation."""
    script_dir = Path(__file__).parent
    base_dir = script_dir.parent
    data_dir = base_dir / "data"
    
    if not data_dir.exists():
        return None
    
    # Find all OJ287_precession_* directories
    dirs = sorted(data_dir.glob("OJ287_precession_*"), reverse=True)
    
    for d in dirs:
        results_file = d / "results.txt"
        if results_file.exists():
            return results_file
    
    return None


def parse_disk_crossings_from_results(filepath):
    """
    Parse disk crossing events from results.txt file.
    
    The results.txt contains crossing times in years from simulation start.
    We map the FIRST crossing to t₀ = JD 2454145.455, then calculate
    all other crossings relative to it.
    
    Parameters:
    -----------
    filepath : str or Path
        Path to results.txt file
    
    Returns:
    --------
    list of dict
        Disk crossing events with JD, year, direction
    """
    crossings = []
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except Exception as e:
        logging.error(f"Error reading results file: {e}")
        return []
    
    # Parse the table: #   Time Before    Time After   Direction    Separation     Orbit #
    # Pattern to match table rows like:
    #    1      6.757200      6.760800          up       70.9671      0.5633
    pattern = r'\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(up|down)\s+([\d.]+)\s+([\d.]+)'
    matches = re.findall(pattern, content)
    
    if not matches:
        logging.warning("No disk crossings found in results file")
        return []
    
    # Parse all crossings
    raw_crossings = []
    for match in matches:
        crossing_num = int(match[0])
        time_before = float(match[1])
        time_after = float(match[2])
        direction = match[3]
        separation = float(match[4])
        orbit_num = float(match[5])
        
        # Average time for the crossing (in years from simulation start)
        avg_time = (time_before + time_after) / 2
        raw_crossings.append({
            'num': crossing_num,
            'time_years': avg_time,
            'direction': direction,
            'separation': separation,
            'orbit': orbit_num
        })
    
    if not raw_crossings:
        return []
    
    # IMPORTANT: Map first crossing to t₀ = JD 2454145.455
    # All other crossings are calculated relative to this
    first_crossing_time = raw_crossings[0]['time_years']
    
    logging.info(f"First crossing in model at t = {first_crossing_time:.3f} years")
    logging.info(f"Mapping to t₀ = JD {T0_JD:.3f}")
    
    for c in raw_crossings:
        # Time difference from first crossing (in years)
        delta_t = c['time_years'] - first_crossing_time
        # Convert to JD relative to t₀
        c['jd'] = T0_JD + delta_t * 365.25
        c['year'] = jd_to_calendar(c['jd'])
        crossings.append(c)
    
    logging.info(f"Parsed {len(crossings)} disk crossings from results file")
    
    return crossings


def calculate_disk_crossings_analytical(num_crossings=30):
    """
    Calculate disk crossing times using analytical approximation
    when results.txt is not available.
    
    Uses typical OJ 287 orbital parameters (~12 year period).
    First crossing is at t₀ = JD 2454145.455.
    
    Parameters:
    -----------
    num_crossings : int
        Number of crossings to generate
    
    Returns:
    --------
    list of dict
        Disk crossing events with JD, year, and direction
    """
    logging.info("Using analytical approximation for disk crossings...")
    logging.info(f"Reference epoch t₀ = JD {T0_JD:.3f}")
    
    # Typical intervals from OJ 287 model (from results.txt pattern)
    # Crossings alternate: ~6.8 years, ~6.4 years, ~7.0 years, etc.
    # Average ~6.7 years between crossings (2 per ~12 year orbit)
    intervals = [6.8, 6.4, 7.0, 6.2, 7.1, 6.4, 6.8, 6.8, 6.8, 6.8, 6.4, 7.0, 6.3, 7.0, 6.4, 6.8]
    
    crossings = []
    
    # First crossing at t₀
    crossings.append({
        'jd': T0_JD,
        'year': jd_to_calendar(T0_JD),
        'direction': 'up',
        'num': 1
    })
    
    # Generate forward crossings
    current_jd = T0_JD
    direction = 'up'
    for i in range(num_crossings):
        interval = intervals[i % len(intervals)] * 365.25  # Convert years to days
        current_jd += interval
        direction = 'down' if direction == 'up' else 'up'
        crossings.append({
            'jd': current_jd,
            'year': jd_to_calendar(current_jd),
            'direction': direction,
            'num': len(crossings) + 1
        })
    
    # Generate backward crossings
    current_jd = T0_JD
    direction = 'up'
    for i in range(num_crossings):
        interval = intervals[i % len(intervals)] * 365.25
        current_jd -= interval
        direction = 'down' if direction == 'up' else 'up'
        crossings.insert(0, {
            'jd': current_jd,
            'year': jd_to_calendar(current_jd),
            'direction': direction,
            'num': 0
        })
    
    # Re-number crossings
    for i, c in enumerate(crossings):
        c['num'] = i + 1
    
    return crossings


def read_photometric_data(filepath):
    """
    Read photometric data from TXT file.
    Format: JD, Magnitude1, Magnitude2
    Supports both European (comma) and American (dot) decimal separators.
    Handles missing values (empty or "-") correctly.
    
    Parameters:
    -----------
    filepath : str
        Path to data file
    
    Returns:
    --------
    tuple
        (jd, mag1, mag2) numpy arrays, or (None, None, None) on error
        mag2 will have NaN where values are missing
    """
    jd = []
    mag1 = []
    mag2 = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Skip header lines (first line or lines starting with #)
        start_idx = 0
        for i, line in enumerate(lines):
            line = line.strip()
            if line and not line.startswith('#'):
                # Check if this line looks like data (starts with a number)
                parts = [p for p in line.split() if p]
                if len(parts) >= 2:
                    try:
                        float(parts[0].replace(',', '.'))
                        # Check if second part is a number (not header text)
                        test_val = parts[1].replace(',', '.')
                        if test_val != '-':
                            float(test_val)
                        start_idx = i
                        break
                    except ValueError:
                        continue
        
        for line in lines[start_idx:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = [p for p in line.split() if p]
            if len(parts) >= 2:
                try:
                    # Convert European format (comma as decimal separator)
                    jd_val = float(parts[0].replace(',', '.'))
                    
                    # Parse second column (mag1)
                    mag1_str = parts[1].replace(',', '.')
                    if mag1_str == '-' or mag1_str == '':
                        mag1_val = np.nan
                    else:
                        mag1_val = float(mag1_str)
                    
                    # Parse third column (mag2) if exists
                    if len(parts) >= 3:
                        mag2_str = parts[2].replace(',', '.')
                        if mag2_str == '-' or mag2_str == '':
                            mag2_val = np.nan
                        else:
                            mag2_val = float(mag2_str)
                    else:
                        mag2_val = np.nan  # No third column
                    
                    jd.append(jd_val)
                    mag1.append(mag1_val)
                    mag2.append(mag2_val)
                    
                except ValueError:
                    continue
                    
    except FileNotFoundError:
        logging.error(f"File not found: {filepath}")
        return None, None, None
    except Exception as e:
        logging.error(f"Error reading file: {e}")
        return None, None, None
    
    if len(jd) == 0:
        logging.error("No valid data found in file!")
        return None, None, None
    
    # Convert to numpy arrays
    jd = np.array(jd)
    mag1 = np.array(mag1)
    mag2 = np.array(mag2)
    
    # Count valid points
    valid_mag1 = np.sum(~np.isnan(mag1))
    valid_mag2 = np.sum(~np.isnan(mag2))
    logging.info(f"Valid mag1 points (стовп. 2): {valid_mag1}")
    logging.info(f"Valid mag2 points (стовп. 3): {valid_mag2}")
    
    return jd, mag1, mag2


def plot_lightcurve(jd, mag1, mag2, filter_name, crossings, output_file='OJ287_lightcurve.png'):
    """
    Plot the light curve with model predictions.
    
    Parameters:
    -----------
    jd : array
        Julian Dates
    mag1 : array
        First magnitude column
    mag2 : array
        Second magnitude column
    filter_name : str
        Filter name for title
    crossings : list of dict
        Disk crossing events with 'jd' and 'direction'
    output_file : str
        Output filename
    """
    fig, ax = plt.subplots(figsize=(16, 9))
    
    # Set white background
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Sort data by JD
    sort_idx = np.argsort(jd)
    jd_sorted = jd[sort_idx]
    mag1_sorted = mag1[sort_idx]
    mag2_sorted = mag2[sort_idx]
    
    # Plot photometric data with different colors
    # Стовпчик 2 (mag1 в коді) - помаранчева
    valid_mag1 = ~np.isnan(mag1_sorted)
    ax.scatter(jd_sorted[valid_mag1], mag1_sorted[valid_mag1], c='#ff6600', s=25, alpha=0.8, 
               edgecolors='black', linewidths=0.3, label='Стовпчик 2', zorder=3)
    # Стовпчик 3 (mag2 в коді) - ніжно голубий, малюється зверху
    valid_mag2 = ~np.isnan(mag2_sorted)
    ax.scatter(jd_sorted[valid_mag2], mag2_sorted[valid_mag2], c='#87CEEB', s=30, alpha=0.9, 
               edgecolors='black', linewidths=0.4, label='Стовпчик 3', zorder=4)
    
    # Invert y-axis (standard for magnitude plots: brighter = smaller magnitude = higher on plot)
    ax.invert_yaxis()
    
    # Get plot limits
    jd_min, jd_max = jd.min(), jd.max()
    jd_range = jd_max - jd_min
    ax.set_xlim(jd_min - 0.02 * jd_range, jd_max + 0.02 * jd_range)
    
 
    
    # Add secondary x-axis with calendar years
    ax2 = ax.secondary_xaxis('top', functions=(jd_to_calendar, calendar_to_jd))
    ax2.set_xlabel('Рік', fontsize=14, color='black', fontweight='bold')
    ax2.tick_params(axis='x', colors='black', labelsize=11)
    
    # Style spines (black)
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)
    for spine in ax2.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)
    
    # Labels and title
    ax.set_xlabel('Юліанська дата (JD)', fontsize=14, color='black', fontweight='bold')
    ax.set_ylabel(f'Магнітуда ({filter_name})', fontsize=14, color='black', fontweight='bold')
    ax.set_title(f'Крива блиску OJ 287 — Фільтр: {filter_name}\n'
                f'Модельні моменти проходження акреційного диска (t₀ = JD {T0_JD:.3f})', 
                fontsize=16, color='black', fontweight='bold', pad=20)
    
    # Grid
    ax.grid(True, alpha=0.3, color='gray', linestyle='-', linewidth=0.5)
    ax.tick_params(axis='both', colors='black', labelsize=11)
    
    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff6600', 
               markersize=10, label='Магнітуда (стовп. 2)', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#87CEEB', 
               markersize=10, label='Магнітуда (стовп. 3)', linestyle='None'),
        Line2D([0], [0], color='#00ff00', linewidth=3.0, 
               label=f't₀ = JD {T0_JD:.3f}'),
        Line2D([0], [0], color='#00ff00', linewidth=2.0, linestyle='-',
               label='Проходження диска'),
    ]
    legend = ax.legend(handles=legend_elements, loc='upper right', 
                      fontsize=10, facecolor='white', edgecolor='black',
                      labelcolor='black')
    
    fig.tight_layout()
    plt.savefig(output_file, dpi=150, facecolor=fig.get_facecolor(), 
                edgecolor='none', bbox_inches='tight')
    plt.close()
    
    logging.info(f"Light curve saved: {output_file}")
    return output_file


def main():
    """Main function."""
    print("=" * 60)
    print("OJ 287 LIGHT CURVE PLOTTER")
    print("with Precession Model Disk Crossing Predictions")
    print("=" * 60)
    print()
    print(f"Reference epoch t₀ = JD {T0_JD:.3f}")
    print()
    
    # Ask for input file path
    filepath = input("Введіть шлях до файлу з фотометричними даними (TXT): ").strip()
    if not filepath:
        print("Помилка: Не вказано шлях до файлу!")
        sys.exit(1)
    
    # Expand user home directory if needed
    filepath = os.path.expanduser(filepath)
    
    if not os.path.isfile(filepath):
        print(f"Помилка: Файл не знайдено: {filepath}")
        sys.exit(1)
    
    # Ask for filter name
    filter_name = input("Введіть назву фільтра (I, V, R, B тощо): ").strip()
    if not filter_name:
        filter_name = "Unknown"
    
    print()
    logging.info(f"Читання даних з: {filepath}")
    logging.info(f"Фільтр: {filter_name}")
    
    # Read photometric data
    jd, mag1, mag2 = read_photometric_data(filepath)
    
    if jd is None:
        print("Помилка: Не вдалося прочитати фотометричні дані!")
        sys.exit(1)
    
    logging.info(f"Прочитано {len(jd)} точок даних")
    logging.info(f"Діапазон JD: {jd.min():.3f} — {jd.max():.3f}")
    logging.info(f"Діапазон років: {jd_to_calendar(jd.min()):.1f} — {jd_to_calendar(jd.max()):.1f}")
    
    # Calculate disk crossings from precession model
    logging.info("Обчислення моментів проходження диска з моделі прецесії...")
    
    # First try to find and read results.txt from precession simulation
    results_file = find_results_file()
    crossings = []
    
    if results_file:
        logging.info(f"Знайдено файл результатів: {results_file}")
        crossings = parse_disk_crossings_from_results(results_file)
    
    if len(crossings) == 0:
        logging.warning("Файл results.txt не знайдено або порожній. Використовую аналітичну модель.")
        crossings = calculate_disk_crossings_analytical(40)
    
    logging.info(f"Модель передбачає {len(crossings)} проходжень диска")
    
    # Filter crossings within data range (with margin)
    jd_margin = 0.15 * (jd.max() - jd.min())
    crossings_in_range = [c for c in crossings 
                         if jd.min() - jd_margin <= c['jd'] <= jd.max() + jd_margin]
    logging.info(f"Проходження в межах графіка: {len(crossings_in_range)}")
    
    # Generate output filename
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    output_file = f"OJ287_lightcurve_{base_name}_{filter_name}.png"
    
    # Plot
    plot_lightcurve(jd, mag1, mag2, filter_name, crossings, output_file)
    
    # Print crossing times
    print()
    print("=" * 70)
    print("МОМЕНТИ ПРОХОДЖЕННЯ ДИСКА (Модельні передбачення)")
    print("=" * 70)
    print(f"{'JD':<16} {'Рік':<10} {'Δt від t₀ (роки)':<18} {'Напрямок':<10}")
    print("-" * 70)
    
    for c in sorted(crossings, key=lambda x: x['jd']):
        if jd.min() - jd_margin <= c['jd'] <= jd.max() + jd_margin:
            delta_t = (c['jd'] - T0_JD) / 365.25
            direction = c.get('direction', '?')
            marker = " <-- t₀" if abs(c['jd'] - T0_JD) < 1 else ""
            direction_ukr = "вгору" if direction == 'up' else "вниз" if direction == 'down' else "?"
            print(f"{c['jd']:<16.3f} {c['year']:<10.2f} {delta_t:<+18.2f} {direction_ukr:<10}{marker}")
    
    print("=" * 70)
    print(f"\nГотово! Крива блиску збережена у: {output_file}")


if __name__ == '__main__':
    main()
