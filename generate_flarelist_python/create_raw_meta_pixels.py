from types import SimpleNamespace
from typing import Union
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.units import Quantity
from sunpy.coordinates import HeliographicStonyhurst
from sunpy.time import TimeRange

from stixpy.calibration.energy import get_elut
from stixpy.calibration.grid import get_grid_transmission
from stixpy.calibration.livetime import get_livetime_fraction
from stixpy.coordinates.frames import STIXImaging
from stixpy.coordinates.transforms import get_hpc_info
from stixpy.io.readers import read_subc_params
from stixpy.config.instrument import STIX_INSTRUMENT, _get_uv_points_data
from stixpy.calibration.visibility import create_visibility, get_elut_correction, _PIXEL_SLICES




from xrayvision.visibility import Visibilities, VisMeta

from stixpy.product.sources import CompressedPixelData, RawPixelData, SummedCompressedPixelData



@u.quantity_input
def create_raw_meta_pixels(
    pixel_data: Union[RawPixelData, CompressedPixelData, SummedCompressedPixelData],
    time_range: Time,
    energy_range: Quantity["energy"],  # noqa: F821
    flare_location: SkyCoord | None = None,
    pixels: str = "top+bot",
    no_shadowing: bool = False,
    return_raw_counts = False,
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
    flare_location
        The coordinates of flare used to calculate grid transmission
    pixels :
        The set of pixels to use to create the meta pixels.
        Allowed values are 'all', 'big' (default), 'small', 'top', 'bottom'.
    no_shadowing : bool optional
        If set to True turn grid shadowing correction off
    Returns
    -------
    `dict`
        Visibility data and sub-collimator information
    """

    # checks if a time bin fully overlaps, is fully within, starts within, or ends within the specified time range.
    pixel_starts = pixel_data.times - pixel_data.duration / 2
    pixel_ends = pixel_data.times + pixel_data.duration / 2

    time_range_start = Time(time_range[0])
    time_range_end = Time(time_range[1])

    t_mask = (
        (pixel_starts >= time_range_start) & (pixel_ends <= time_range_end)  # fully within the time range
        | (time_range_start <= pixel_starts) & (time_range_end >= pixel_ends)  #  fully overlaps the time range
        | (pixel_starts <= time_range_start) & (pixel_ends >= time_range_start)  # starts within the time range
        | (pixel_starts <= time_range_end) & (pixel_ends >= time_range_end)  #  ends within the time range
    )

    e_mask = (pixel_data.energies["e_low"] >= energy_range[0]) & (pixel_data.energies["e_high"] <= energy_range[1])

    t_ind = np.argwhere(t_mask).ravel()
    e_ind = np.argwhere(e_mask).ravel()

    time_range = TimeRange(
        pixel_data.times[t_ind[0]] - pixel_data.duration[t_ind[0]] / 2,
        pixel_data.times[t_ind[-1]] + pixel_data.duration[t_ind[-1]] / 2,
    )

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

    _, livefrac, _ = get_livetime_fraction(triggers / pixel_data.data["timedel"].to("s").reshape(-1, 1))

    pixel_data.data["livefrac"] = livefrac

    # e_cor_high, e_cor_low = get_elut_correction(e_ind, pixel_data)


    # Get counts and other data.
    idx_pix = _PIXEL_SLICES.get(pixels.lower(), None)
    if idx_pix is None:
        raise ValueError(f"Unrecognised input for 'pixels': {pixels}. Supported values: {list(_PIXEL_SLICES.keys())}")
    counts = pixel_data.data["counts"].astype(float)
    count_errors = np.sqrt(pixel_data.data["counts_comp_err"].astype(float).value ** 2 + counts.value) * u.ct
    ct = counts[t_ind][..., idx_pix, e_ind]
    ct[..., 0] = ct[..., 0] 
    # * e_cor_low[..., idx_pix]
    ct[..., -1] = ct[..., -1]
    # * e_cor_high[..., idx_pix]
    ct_error = count_errors[t_ind][..., idx_pix, e_ind]
    ct_error[..., 0] = ct_error[..., 0] 
    # * e_cor_low[..., idx_pix]
    ct_error[..., -1] = ct_error[..., -1] 
    # * e_cor_high[..., idx_pix]


    lt = (livefrac * pixel_data.data["timedel"].reshape(-1, 1).to("s"))[t_ind].sum(axis=0)

    ct_summed = ct.sum(axis=(0, 3))  # .astype(float)
    ct_error_summed = np.sqrt(np.sum(ct_error**2, axis=(0, 3)))


    # can feactor below to its separate function

    # returns raw counts for each of 24 subcolimators, 4 bottom counts, 4 top counts
    # with no elut correction

    raw_counts = None
    raw_counts_24 = None
    if return_raw_counts and pixels == "top+bot":
        

        # shape is (356,32,8,17)
        print(f'counts shape: {counts.shape}') 


        # shape is (1,32,8,12)
        print(f'ct shape: {ct.shape}')         

        #ct_summed = ct.sum(axis=(0, 3))  # sums over tiem ranges and energy ranges
        #(32,8)
        print(f'ct_summed shape: {ct_summed.shape}') 

        abcd_top_counts = ct_summed[:, :4]
        # (32, 4) 
        abcd_bot_counts = ct_summed[:, 4:]

       
        print(f'abcd_top_counts shape: {abcd_top_counts.shape}')

        print(f'abcd_top_counts example collimator 0: {abcd_top_counts[:1]}')
        print(f'abcd_bot_counts shape: {abcd_bot_counts.shape}')

        print(f'abcd_bot_counts example collimator 0: {abcd_bot_counts[:1]}')


        # stack them together to get shape 2,32,4
        raw_counts = np.stack([abcd_top_counts, abcd_bot_counts], axis=0) 
        print(f'raw counts shape: {raw_counts.shape}')


       

        # better to do this step outside of the function i think when getting meta pixels
        # raw_counts_24 = raw_counts[:, isc_24, :]
        #                   take all, take isc 24, take all

        # print(f'final abcd raw counts 24 shape: {raw_counts_24.shape}')
    

    ## im returning raw counts fix above
      

    

















        


    if not no_shadowing:
        if flare_location is None or not isinstance(flare_location, SkyCoord):
            raise ValueError("flare_location must be a SkyCoord object if using grid shadowing correction.")
        if not isinstance(flare_location.frame, STIXImaging) and flare_location.obstime != time_range.center:
            roll, solo_heeq, stix_pointing = get_hpc_info(time_range.start, time_range.end)
            flare_location = flare_location.transform_to(STIXImaging(obstime=time_range.center, observer=solo_heeq))
        grid_shadowing = get_grid_transmission(flare_location)
        ct_summed = ct_summed / grid_shadowing.reshape(-1, 1) / 4  # transmission grid ~ 0.5*0.5 = .25
        ct_error_summed = ct_error_summed / grid_shadowing.reshape(-1, 1) / 4

    abcd_counts = ct_summed.reshape(ct_summed.shape[0], -1, 4).sum(axis=1)
    abcd_count_errors = np.sqrt((ct_error_summed.reshape(ct_error_summed.shape[0], -1, 4) ** 2).sum(axis=1))

    abcd_rate = abcd_counts / lt.reshape(-1, 1)
    abcd_rate_error = abcd_count_errors / lt.reshape(-1, 1)

    e_bin = pixel_data.energies[e_ind][-1]["e_high"] - pixel_data.energies[e_ind][0]["e_low"]
    abcd_rate_kev = abcd_rate / e_bin
    abcd_rate_error_kev = abcd_rate_error / e_bin

    pixel_areas = STIX_INSTRUMENT.pixel_config["Area"].to("cm2")
    areas = pixel_areas[idx_pix].reshape(-1, 4).sum(axis=0)
    

    meta_pixels = {
        "abcd_rate_kev": abcd_rate_kev,
        "abcd_rate_kev_cm": abcd_rate_kev / areas,
        "abcd_rate_error_kev_cm": abcd_rate_error_kev / areas,
        "time_range": time_range,
        "energy_range": energy_range,
        "pixels": pixels,
        "areas": areas,
    }

    if return_raw_counts:
        print(" im returning")
        return meta_pixels, raw_counts # or raw counts hmm
    
    return meta_pixels


