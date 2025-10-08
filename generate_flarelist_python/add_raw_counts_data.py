from typing import Union


from astropy.units import Quantity


import pandas as pd
import numpy as np
from astropy.time import Time
from astropy import units as u
from sunpy.time import parse_time
from stixpy.product import Product


import logging

from stixpy.calibration.livetime import get_livetime_fraction

from stixpy.config.instrument import STIX_INSTRUMENT
from stixpy.product.sources import CompressedPixelData, RawPixelData, SummedCompressedPixelData

def add_raw_counts_data(flare_list_with_files, save_csv=False):
    """
        Adds 24 sub-collimators raw counts, top and bottom, a to d, as separate column

        Parameters
        ----------
        flare_list_with_files : pd.DataFrame
            DataFrame containing flare information including file paths (`filenames`) to associated `.fits` files.

        """

    logging.info('Getting raw counts and adding them...')

    results = {}

    # for now hardcoded, can be refactored to be passed as argument
    # we only care about the 24 sub-collimators below
    # -1 for 0-index
    isc_24 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]) - 1

    # create all raw count columns
    for sc in isc_24:
        for letter in ["a", "b", "c", "d"]:
            results[f'{sc + 1}_{letter}_top'] = []
            results[f'{sc + 1}_{letter}_bot'] = []

    results['attenuator'] = []

    for i, row in flare_list_with_files.iterrows():
        energy_range = [4, 10] * u.keV

        # Define a 20s time range around peak time
        tstart = parse_time(row["peak_UTC"]) - 20 * u.s
        tend = parse_time(row["peak_UTC"]) + 20 * u.s
        time_range = [tstart.strftime("%Y-%m-%dT%H:%M:%S"), tend.strftime("%Y-%m-%dT%H:%M:%S")]
        cpd_file = row["filenames"]
        att = False  # Default value for attenuator


        cpd_sci = Product(cpd_file)

        # Check for attenuator status by looking for any 'rcr' data points in the time range
        # as the att_in column in the operational flarelist isnt working.
        if np.any(cpd_sci.data[(cpd_sci.data["time"] >= tstart) & (cpd_sci.data["time"] <= tend)]["rcr"]):
            att = True
            energy_range = [12, 25] * u.keV
        # print(att)

        try:
            raw_counts = get_raw_counts(cpd_sci, time_range, energy_range)

            # filter to only the 24 sub-collimators we care about
            raw_counts_filtered = raw_counts[:, isc_24, :]

            for sc_idx, sc in enumerate(isc_24):
                for abcd_idx, letter in enumerate(["a", "b", "c", "d"]):
                    results[f'{sc + 1}_{letter}_top'].append(raw_counts_filtered[0, sc_idx, abcd_idx].value)
                    results[f'{sc + 1}_{letter}_bot'].append(raw_counts_filtered[1, sc_idx, abcd_idx].value)

        except Exception as e:
            logging.error(f'error getting raw counts for flare {i}: {e}')
            # Append NaN for all raw counts
            for sc in isc_24:
                for letter in ["a", "b", "c", "d"]:
                    results[f'{sc + 1}_{letter}_top'].append(np.nan)
                    results[f'{sc + 1}_{letter}_bot'].append(np.nan)

        results["attenuator"].append(att)

    results = pd.DataFrame(results)
    flare_list_with_raw_counts = pd.concat([flare_list_with_files.reset_index(drop=True), results], axis=1)

    times_flares = pd.to_datetime(flare_list_with_raw_counts["peak_UTC"])

    if save_csv:
        filename = f"stix_flarelist_w_raw_counts_{times_flares.min().strftime('%Y%m%d')}_{times_flares.max().strftime('%Y%m%d')}.csv"
        flare_list_with_raw_counts.to_csv(filename, index=False, index_label=False)
        logging.info(f'Saved flare list to {filename}')

    return flare_list_with_raw_counts

@u.quantity_input
def get_raw_counts(
        pixel_data: Union[RawPixelData, CompressedPixelData, SummedCompressedPixelData],
        time_range: Time,
        energy_range: Quantity["energy"],  # noqa: F821
):
    r"""
    Create meta-pixels by summing data with in given time and energy range.

    Parameters
    ----------
    pixel_data
        Input pixel data
    time_range :
        Start and end times
    energy_range
        Start and end energies
    Returns
    -------
    `np array`
        raw sub-collimator counts top and bottom shape (2,32,4)
    """

    # checks if a time bin fully overlaps, is fully within, starts within, or ends within the specified time range.
    pixel_starts = pixel_data.times - pixel_data.duration / 2
    pixel_ends = pixel_data.times + pixel_data.duration / 2

    time_range_start = Time(time_range[0])
    time_range_end = Time(time_range[1])

    t_mask = (
            (pixel_starts >= time_range_start) & (pixel_ends <= time_range_end)  # fully within the time range
            | (time_range_start <= pixel_starts) & (time_range_end >= pixel_ends)  # fully overlaps the time range
            | (pixel_starts <= time_range_start) & (pixel_ends >= time_range_start)  # starts within the time range
            | (pixel_starts <= time_range_end) & (pixel_ends >= time_range_end)  # ends within the time range
    )

    e_mask = (pixel_data.energies["e_low"] >= energy_range[0]) & (pixel_data.energies["e_high"] <= energy_range[1])

    t_ind = np.argwhere(t_mask).ravel()
    e_ind = np.argwhere(e_mask).ravel()


    changed = []
    for column in ["rcr", "pixel_masks", "detector_masks"]:
        if np.unique(pixel_data.data[column][t_ind], axis=0).shape[0] != 1:
            changed.append(column)
    if len(changed) > 0:
        raise ValueError(
            f"The following: {', '.join(changed)} changed in the selected time interval "
            f"please select a time interval where these are constant."
        )

    trigger_to_detector = STIX_INSTRUMENT.subcol_adc_mapping

    # Map the triggers to all 32 detectors
    triggers = pixel_data.data["triggers"][:, trigger_to_detector].astype(float)[...]

    livefrac, *_ = get_livetime_fraction(triggers / pixel_data.data["timedel"].to("s").reshape(-1, 1))

    pixel_data.data["livefrac"] = livefrac


    # get top and bottom
    idx_pix = slice(0, 8)
    counts = pixel_data.data["counts"].astype(float)
    ct = counts[t_ind][..., idx_pix, e_ind]



    ct_summed = ct.sum(axis=(0, 3))  # .astype(float)

    # returns raw counts for each of 32 sub-collimators, 4 bottom counts, 4 top counts
    raw_counts = None


    # (32, 4)
    abcd_top_counts = ct_summed[:, :4]
    abcd_bot_counts = ct_summed[:, 4:]

    # stack them together to get shape 2,32,4
    raw_counts = np.stack([abcd_top_counts, abcd_bot_counts], axis=0)




    return raw_counts
