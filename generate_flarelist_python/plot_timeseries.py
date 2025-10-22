import numpy as np
import pandas as pd
from astropy.time import Time
from stixpy.product import Product
import logging
import matplotlib.pyplot as plt
from matplotlib.dates import AutoDateLocator, ConciseDateFormatter
from pathlib import Path
from datetime import datetime
import os

from astropy import units as u



def do_plotting(flare_list_with_files, output_dir='flare_plots'):
    """
    Plot flares grouped by day with all flares from each day in one figure.
    
    Parameters
    ----------
    flare_list_with_files : pd.DataFrame
        DataFrame with flare raw_fits including 'filenames', 'peak_UTC', 'flare_id', 'att_in'
    output_dir : str
        Directory to save the plots
    """
    
    logging.info("Creating plots for flares grouped by day")
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Extract date from peak_UTC and group by day
    flare_list_with_files['date'] = pd.to_datetime(flare_list_with_files['peak_UTC']).dt.date
    grouped_by_day = flare_list_with_files.groupby('date')
    
    print(f"Processing {len(grouped_by_day)} days with flares")
    
    for date, day_flares in grouped_by_day:
        print(f"\n{'='*60}")
        print(f"Processing date: {date} ({len(day_flares)} flares)")
        print(f"{'='*60}")
        
        # Calculate grid dimensions for subplots
        n_flares = len(day_flares)
        n_cols = min(3, n_flares)  # Max 3 columns
        n_rows = int(np.ceil(n_flares / n_cols))
        
        # Create figure with subplots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(8*n_cols, 5*n_rows))
        if n_flares == 1:
            axes = np.array([axes])  # Make it iterable
        axes = axes.flatten()
        
        flare_count = 0
        for idx, (_, row) in enumerate(day_flares.iterrows()):
            ax = axes[idx]
            
            cpd_file = row["filenames"]
            flare_id = row["flare_id"]
            
            # Skip if file has issues
            if pd.isna(cpd_file) or cpd_file == '' or 'file_issue' in str(cpd_file).lower():
                print(f"  Skipping flare {flare_id}: file issue")
                ax.text(0.5, 0.5, f'Flare {flare_id}\nFile Issue', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            
            try:
                # Load product
                cpd_sci = Product(cpd_file)
                
                # Get flare peak time
                peak_time = Time(row["peak_UTC"])
                
                # Use att_in from CSV or check from raw_fits
                att = bool(row.get("att_in", False))
                
                # Find ALL time periods where attenuator is IN across the entire FITS file
                rcr_mask = cpd_sci.data["rcr"] > 0
                att_times = cpd_sci.data[rcr_mask]["time"]
                
                # Set energy range based on attenuator status
                if att:
                    energy_range = [12, 25] * u.keV
                    plot_color = 'orange'
                else:
                    energy_range = [4, 10] * u.keV
                    plot_color = 'blue'
                
                e_mask = (cpd_sci.energies["e_low"] >= energy_range[0]) & (cpd_sci.energies["e_high"] <= energy_range[1])
                e_ind = np.argwhere(e_mask).ravel()
                
                if len(e_ind) == 0:
                    print(f"  Skipping flare {flare_id}: no raw_fits in energy range")
                    ax.text(0.5, 0.5, f'Flare {flare_id}\nNo raw_fits in energy range',
                           ha='center', va='center', transform=ax.transAxes, fontsize=10)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue
                
                print(f"  Plotting flare {flare_id} (att={'IN' if att else 'OUT'})")
                
                # Plot timeseries
                e_ind_summed = [[e_ind[0], e_ind[-1]]]
                cpd_sci.plot_timeseries(energy_indices=e_ind_summed, axes=ax, color=plot_color)
                
                # Add vertical line for flare peak
                ax.axvline(peak_time.to_datetime(), color='red', linestyle='--', linewidth=2,
                          label='Peak', zorder=5)
                
                # Shade attenuator periods
                ylim = ax.get_ylim()
                if len(att_times) > 0:
                    att_indices = np.where(rcr_mask)[0]
                    if len(att_indices) > 0:
                        splits = np.where(np.diff(att_indices) > 1)[0] + 1
                        att_groups = np.split(att_indices, splits)
                        
                        first_group = True
                        for group in att_groups:
                            if len(group) > 0:
                                start_time = cpd_sci.times[group[0]].to_datetime()
                                end_time = cpd_sci.times[group[-1]].to_datetime()
                                if first_group:
                                    ax.axvspan(start_time, end_time, alpha=0.15, color='orange', 
                                             zorder=0, label='ATT IN')
                                    first_group = False
                                else:
                                    ax.axvspan(start_time, end_time, alpha=0.15, color='orange', zorder=0)
                ax.set_ylim(ylim)
                
                # Set title
                att_status = "ATT IN" if att else "ATT OUT"
                ax.set_title(f'Flare {flare_id} | {att_status}\n{energy_range[0].value}-{energy_range[1].value} keV',
                           fontsize=10)
                ax.legend(loc='upper right', fontsize=8)
                
                # Set up date formatting
                locator = AutoDateLocator()
                formatter = ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
                
                flare_count += 1
                
            except Exception as e:
                print(f"  Error plotting flare {flare_id}: {e}")
                ax.text(0.5, 0.5, f'Flare {flare_id}\nError: {str(e)[:50]}...', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=9, wrap=True)
                ax.set_xticks([])
                ax.set_yticks([])
        
        # Hide empty subplots
        for idx in range(n_flares, len(axes)):
            axes[idx].set_visible(False)
        
        # Set overall title
        fig.suptitle(f'Flares on {date} ({flare_count}/{n_flares} plotted successfully)', 
                    fontsize=16, fontweight='bold')
        fig.tight_layout()
        
        # Save figure
        output_file = Path(output_dir) / f'flares_{date}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
        
        plt.close(fig)
    
    print(f"\n{'='*60}")
    print(f"All plots saved to: {output_dir}/")
    print(f"{'='*60}")






if __name__ == "__main__":
    print("Reading flare raw_fits...")
    flare_list_with_files = pd.read_csv('stix_flarelist_w_locations_20240901_20240909.csv')
    print(f"Loaded {len(flare_list_with_files)} flares")
    
    # Create plots grouped by day
    do_plotting(flare_list_with_files, output_dir='flare_plots')





