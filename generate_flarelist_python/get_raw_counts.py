from typing import Union


import astropy.units as u
import numpy as np
from astropy.time import Time
from astropy.units import Quantity


from stixpy.calibration.livetime import get_livetime_fraction

from stixpy.config.instrument import STIX_INSTRUMENT
from stixpy.product.sources import CompressedPixelData, RawPixelData, SummedCompressedPixelData

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
