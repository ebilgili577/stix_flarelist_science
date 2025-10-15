import numpy as np
import pandas as pd
from astropy.time import Time
from stixpy.product import Product
import logging
import matplotlib.pyplot as plt
from matplotlib.dates import AutoDateLocator, ConciseDateFormatter


from astropy import units as u



def do_plotting(flare_list_with_files, start_index=0):

    logging.info("creating plots for flares")

    # skip to flare x
    flare_items = list(flare_list_with_files.iterrows())
    if start_index > 0:
        print(f"Skipping to flare #{start_index} (out of {len(flare_items)} total flares)")
        flare_items = flare_items[start_index:]

    for i, row in flare_items:

        cpd_file = row["filenames"]
        att = False  # Default value for attenuator




        cpd_sci = Product(cpd_file)

        # get flare peak time
        peak_time = Time(row["peak_UTC"])

        # Check for attenuator status around the flare peak
        time_window = 20 * u.s  # check +/- 20 sec around peak
        tstart_check = peak_time - time_window
        tend_check = peak_time + time_window
        time_mask_att = (cpd_sci.data["time"] >= tstart_check) & (cpd_sci.data["time"] <= tend_check)
        if np.any(cpd_sci.data[time_mask_att]["rcr"]):
            att = True

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

        print(f'Plotting flare: {row["flare_id"]}')

        # Create single interactive plot with full FITS range
        fig = plt.figure(figsize=(14, 7))
        ax = fig.add_subplot(111)

        # plotting using cpd timeseries

        # sum over energy ranges 4-10
        e_ind_summed = [[e_ind[0], e_ind[-1]]]
        cpd_sci.plot_timeseries(energy_indices=e_ind_summed, axes=ax, color=plot_color)

        # vertical line for flares peak
        ax.axvline(peak_time.to_datetime(), color='red', linestyle='--', linewidth=2,
                   label=f'Flare Peak', zorder=5)




        # bg shadin for att in
        ylim = ax.get_ylim()
        if len(att_times) > 0:
            # Get start and end times for each data point where rcr > 0
            att_indices = np.where(rcr_mask)[0]
            if len(att_indices) > 0:
                # Group consecutive indices to find contiguous attenuator periods
                splits = np.where(np.diff(att_indices) > 1)[0] + 1
                att_groups = np.split(att_indices, splits)

                # Shade each contiguous period
                first_group = True
                for group in att_groups:
                    if len(group) > 0:
                        start_time = cpd_sci.times[group[0]].to_datetime()
                        end_time = cpd_sci.times[group[-1]].to_datetime()
                        if first_group:
                            ax.axvspan(start_time, end_time, alpha=0.15, color='orange', zorder=0, label='Attenuator IN')
                            first_group = False
                        else:
                            ax.axvspan(start_time, end_time, alpha=0.15, color='orange', zorder=0)
        ax.set_ylim(ylim)  # Restore y-limits

        att_status = "ATT IN" if att else "ATT OUT"
        ax.set_title(f'Full FITS Range - Flare: {row["flare_id"]} | {att_status}\n'
                     f'Energy: {energy_range[0].value}-{energy_range[1].value} keV | Use toolbar to zoom/pan')
        ax.legend(loc='upper right')

        # Set up automatic date formatting that adapts to zoom level
        locator = AutoDateLocator()
        formatter = ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)

        fig.tight_layout()

        plt.show()
        plt.close('all')






if __name__ == "__main__":
    print("reading data")
    flare_list_with_files = pd.read_csv('test_flare_data.csv')
    do_plotting(flare_list_with_files)





