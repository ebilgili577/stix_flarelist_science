# Solar Orbiter/STIX science flare list 

This is a repository for a study of the STIX flare list using the scientific pixel data.
This flarelist builds upon the operational STIX flare list that is avalable through [stixdcpy API](https://github.com/i4Ds/stixdcpy).
From the events in this list with available pixel data, and with counts above 1000 in the 4-10 keV energy band an image was generated and the location of the flare estimated. 
Currently it consists of ~25150 flares over the time period of 2021-01-01 to 2025-02-28.

 **just FYI I need to update this** 
 
Here we provide the flare list with the coordinates of the flare estimated, and in several coordinate frames and information whether that flare was observed from Earth.
The flarelist if provided in this file: `STIX_flarelist_w_locations_20210318_20250228_version1_python.csv`

This can be read in python using `pd.read_csv`
e.g. 

```
>>> import pandas as pd
>>> stix_flarelist = pd.read_csv("STIX_flarelist_w_locations_20210318_20250228_version1_python.csv")`
```
or similary in IDL using `READ_CSV()`.

## STIX flare list with raw counts:

In this file, the flarelist contains 213 columns with the following information:

### Basic Flare Information
* `flare_id` : Unique ID of flare from stixdcpy
* `start_UTC` : Start time of flare (UTC)
* `end_UTC` : End time of flare (UTC)
* `peak_UTC` : Peak time of flare (UTC)
* `duration` : Duration of flare in seconds
* `att_in` : Boolean indicating whether the attenuator was in or not
  - **Important:** Energy range selection depends on attenuator state:
    - `att_in = False` (normal): 4-10 keV energy range used for analysis
    - `att_in = True` (attenuator in): 12-25 keV energy range used for analysis

### Raw Pixel Counts (192 columns)
Raw counts for each detector-subcollimator combination:
* `{detector}_{subcollimator}_top` : Raw counts for top detector
* `{detector}_{subcollimator}_bottom` : Raw counts for bottom detector
  - **Detectors:** 0-7, 16-31 (32 detectors total)
  - **Subcollimators:** a, b, c, d (4 per detector)
  - **Example columns:** `0_a_top`, `0_a_bottom`, `1_b_top`, `16_c_bottom`, `31_d_top`, etc.
  - **Total:** 32 detectors × 4 subcollimators × 2 (top/bottom) = 256 potential columns
    - Note: Only detectors 0-7 and 16-31 are present (detectors 8-15 excluded)

### Location Information
* `loc_x_stix` : X location in STIX pixel coordinates
* `loc_y_stix` : Y location in STIX pixel coordinates
* `CFL_X_arcsec` : Coarse Flare Locator X position in arcsec
* `CFL_Y_arcsec` : Coarse Flare Locator Y position in arcsec
* `sidelobes_ratio` : Ratio of sidelobes in imaging analysis
* `error` : Error estimate from location fitting

### GOES Data
* `GOES_class` : GOES class at time of flare (e.g., C1.5, M2.0, X1.0)
* `GOES_flux` : GOES 1-8 Å flux at time of flare (W/m²)
* `goes_estimated_min_class` : Estimated minimum GOES class based on STIX data
* `goes_estimated_max_class` : Estimated maximum GOES class based on STIX data
* `goes_estimated_mean_class` : Estimated mean GOES class based on STIX data
* `goes_estimated_min_flux` : Estimated minimum GOES flux based on STIX data
* `goes_estimated_max_flux` : Estimated maximum GOES flux based on STIX data
* `goes_estimated_mean_flux` : Estimated mean GOES flux based on STIX data

### Additional Metadata
* `filenames` : Associated STIX data filenames used for analysis 


How it was generated:
--------------------
Overview

This pipeline processes STIX flare observations from the Solar Orbiter mission, providing a complete end-to-end analysis framework to:

* Retrieve operational flare lists from the STIX Data Center.
* Filter and associate files based on intensity thresholds.
* Estimate flare locations and attenuator status using STIX imaging.
* Process positional calculations and transform coordinates to different frames.
* Check visibility of flares from Earth and save the final processed flare list to a CSV file.


Setup of dependencies
--------------------

This project requires Python 3.12 and various scientific libraries (astropy, sunpy, stixpy, etc.).

### Installation with Conda

1. **Create the environment:**
   ```bash
   conda env create -f environment.yml
   ```

2. **Activate the environment:**
   ```bash
   conda activate stix_flarelist
   ```

3. **Run scripts:**
   ```bash
   cd generate_flarelist_python
   python flarelist_generate.py
   ```

### Available Scripts

- `get_operational_flarelist.py` - Retrieve operational flare list from STIX
- `stx_estimate_flare_location.py` - Estimate flare locations
- `add_raw_counts_data.py` - Add raw count data to flare list
- `plot_timeseries.py` - Visualize time series



### Update Environment

If new dependencies are added:

```bash
conda env update -f environment.yml --prune
``` 
