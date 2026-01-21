#!/usr/bin/env python3
"""
OJ287 Light Curve Visualization with Disk Crossing Events
Plots observational data and marks predicted disk crossing times
from the CBwaves binary black hole model.
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path

# === Configuration ===
# Reference epoch: JD at model time = 0
# Calibrated so that disk crossings align with known OJ287 outbursts
# Based on Valtonen et al. papers: 2007 Sep outburst corresponds to crossing #16
# Adjusted to match observed outburst at JD ~2454350 (Sep 2007)
EPOCH_JD = 2415207.5  # Calibrated for OJ287

# Additional known OJ287 outburst predictions from literature
# These are added manually since the model time range may not cover them all
ADDITIONAL_OUTBURSTS = [
    {'year': 2019.57, 'jd': 2458666, 'label': 'July 2019', 'source': 'Valtonen+2019'},
]

# === Parse observational data ===
def parse_observations(filepath):
    """Parse JD and magnitude from observation file (European number format)."""
    jd_list = []
    mag_list = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    for line in lines[2:]:  # Skip header lines
        line = line.strip()
        if not line:
            continue
        
        # Split by whitespace
        parts = line.split()
        if len(parts) >= 2:
            try:
                # Convert European format (comma as decimal separator)
                jd = float(parts[0].replace(',', '.'))
                mag = float(parts[1].replace(',', '.'))
                jd_list.append(jd)
                mag_list.append(mag)
            except ValueError:
                # Skip lines with missing or invalid data
                continue
        elif len(parts) == 1:
            # Some lines might have only JD without mag
            continue
    
    return np.array(jd_list), np.array(mag_list)


def parse_disk_crossings(filepath):
    """Parse disk crossing events from results.txt file."""
    crossings = []
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Find the table section
    pattern = r'\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(up|down)\s+([\d.]+)\s+([\d.]+)'
    matches = re.findall(pattern, content)
    
    for match in matches:
        crossing_num = int(match[0])
        time_before = float(match[1])
        time_after = float(match[2])
        direction = match[3]
        separation = float(match[4])
        orbit_num = float(match[5])
        
        # Average time for the crossing
        avg_time = (time_before + time_after) / 2
        crossings.append({
            'num': crossing_num,
            'time_years': avg_time,
            'direction': direction,
            'separation': separation,
            'orbit': orbit_num
        })
    
    return crossings


def model_time_to_jd(time_years, epoch_jd=EPOCH_JD):
    """Convert model time (years) to Julian Date."""
    return epoch_jd + time_years * 365.25


def jd_to_calendar(jd):
    """Convert JD to approximate calendar year."""
    return 2000.0 + (jd - 2451545.0) / 365.25


def find_latest_data_dir(base_dir):
    """Find the most recent OJ287_precession_* directory."""
    data_dir = base_dir / "data"
    dirs = sorted(data_dir.glob("OJ287_precession_*"), reverse=True)
    if not dirs:
        raise FileNotFoundError(f"No OJ287_precession_* directories found in {data_dir}")
    return dirs[0]


def main():
    # File paths
    script_dir = Path(__file__).parent
    base_dir = script_dir.parent
    
    # Auto-detect latest data directory
    data_subdir = find_latest_data_dir(base_dir)
    print(f"Using data directory: {data_subdir.name}")
    
    obs_file = Path("/Users/kinnamou/Desktop/I фільтр .txt")
    results_file = data_subdir / "results.txt"
    
    # Parse data
    print("Reading observational data...")
    jd_obs, mag_obs = parse_observations(obs_file)
    print(f"  Loaded {len(jd_obs)} observations")
    print(f"  JD range: {jd_obs.min():.1f} - {jd_obs.max():.1f}")
    print(f"  Year range: {jd_to_calendar(jd_obs.min()):.1f} - {jd_to_calendar(jd_obs.max()):.1f}")
    print(f"  Magnitude range: {mag_obs.min():.2f} - {mag_obs.max():.2f}")
    
    print("\nReading disk crossing data...")
    crossings = parse_disk_crossings(results_file)
    print(f"  Found {len(crossings)} disk crossing events")
    
    # Convert crossing times to JD
    for c in crossings:
        c['jd'] = model_time_to_jd(c['time_years'])
        c['year'] = jd_to_calendar(c['jd'])
    
    # Filter crossings to observation range (with larger margin to show more events)
    jd_min, jd_max = jd_obs.min() - 1000, jd_obs.max() + 500
    visible_crossings = [c for c in crossings if jd_min <= c['jd'] <= jd_max]
    
    print(f"\nDisk crossings in observation window:")
    for c in visible_crossings:
        print(f"  #{c['num']}: JD {c['jd']:.1f} ({c['year']:.2f}) - {c['direction']}")
    
    # === Create the plot ===
    fig, ax = plt.subplots(figsize=(16, 9))
    
    # Set dark background for dramatic effect
    fig.patch.set_facecolor('#0a0a1a')
    ax.set_facecolor('#0a0a1a')
    
    # Sort observations by JD for plotting
    sort_idx = np.argsort(jd_obs)
    jd_sorted = jd_obs[sort_idx]
    mag_sorted = mag_obs[sort_idx]
    
    # Convert JD to years for x-axis
    years_obs = jd_to_calendar(jd_sorted)
    
    # Plot disk crossing events FIRST (behind data)
    # Get y-axis limits based on data
    mag_min, mag_max = mag_sorted.min() - 0.3, mag_sorted.max() + 0.3
    
    # Single color for all crossings
    crossing_color = '#00ffcc'
    
    for c in visible_crossings:
        year = c['year']
        
        # Draw vertical line for disk crossing
        ax.axvline(x=year, color=crossing_color, alpha=0.7, linestyle='-', linewidth=2.5, zorder=1)
    
    # Plot observations as scatter points
    scatter = ax.scatter(years_obs, mag_sorted, 
                        c=mag_sorted, cmap='plasma_r',
                        s=35, alpha=0.85, edgecolors='white', linewidths=0.4,
                        zorder=3, label='Спостереження (I фільтр)')
    
    # Plot additional known outbursts from literature (same style as model crossings)
    for outburst in ADDITIONAL_OUTBURSTS:
        year = outburst['year']
        ax.axvline(x=year, color=crossing_color, alpha=0.7, linestyle='-', linewidth=2.5, zorder=1)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('Магнітуда', fontsize=12, color='white')
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
    
    # Style the plot
    ax.set_xlabel('Рік', fontsize=14, color='white', fontweight='bold')
    ax.set_ylabel('Магнітуда (I фільтр)', fontsize=14, color='white', fontweight='bold')
    ax.set_title('Крива блиску OJ 287 та моменти проходження акреційного диску\n'
                'Binary Black Hole Model - CBwaves', 
                fontsize=16, color='white', fontweight='bold', pad=20)
    
    # Set y-axis limits (inverted: higher values at bottom = fainter)
    ax.set_ylim(mag_max + 0.2, mag_min - 0.2)
    
    # Add secondary x-axis for JD
    ax2 = ax.secondary_xaxis('top', functions=(
        lambda x: 2451545.0 + (x - 2000.0) * 365.25,  # year to JD
        lambda jd: 2000.0 + (jd - 2451545.0) / 365.25  # JD to year
    ))
    ax2.set_xlabel('Юліанська дата (JD)', fontsize=12, color='white', fontweight='bold')
    ax2.tick_params(axis='x', colors='white', labelsize=10)
    
    # Grid
    ax.grid(True, alpha=0.2, color='white', linestyle='-', linewidth=0.5)
    
    # Spine colors
    for spine in ax.spines.values():
        spine.set_color('white')
        spine.set_linewidth(1.5)
    for spine in ax2.spines.values():
        spine.set_color('white')
        spine.set_linewidth(1.5)
    
    # Tick colors
    ax.tick_params(axis='both', colors='white', labelsize=11)
    
    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff9500', 
               markersize=10, label='Спостереження (I фільтр)', linestyle='None'),
        Line2D([0], [0], color='#00ffcc', linewidth=2.5, alpha=0.7,
               label='Проходження диску'),
    ]
    legend = ax.legend(handles=legend_elements, loc='lower right', 
                      fontsize=10, facecolor='#1a1a2e', edgecolor='white',
                      labelcolor='white')
    
    # Add info text
    info_text = f"Епоха моделі: JD {EPOCH_JD:.1f}\nОрбітальний період: ~12 років"
    ax.text(0.02, 0.02, info_text, transform=ax.transAxes,
            fontsize=9, color='#888888', va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#0a0a1a', 
                     edgecolor='#444444', alpha=0.9))
    
    plt.tight_layout()
    
    # Save figure
    output_path = data_subdir / "OJ287_lightcurve_with_crossings.png"
    plt.savefig(output_path, dpi=150, facecolor=fig.get_facecolor(), edgecolor='none',
                bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")
    
    plt.close(fig)


if __name__ == '__main__':
    main()
