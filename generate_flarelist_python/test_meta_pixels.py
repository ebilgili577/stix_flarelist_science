import warnings
from sunpy.util import SunpyDeprecationWarning
from add_raw_counts_data import add_raw_counts_data

# TODO: Check if this is needed
warnings.filterwarnings("ignore", category=SunpyDeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sunpy")
warnings.filterwarnings("ignore", category=UserWarning, module="stixpy")
warnings.filterwarnings("ignore", category=UserWarning, module="astropy")

# Suppress stixpy logging completely
import logging
logging.getLogger('stixpy').setLevel(logging.CRITICAL)
logging.getLogger('stixpy.coordinates.transforms').setLevel(logging.CRITICAL)

from flarelist_generate import fetch_operational_flare_list, filter_and_associate_files
from astropy.time import Time
from tqdm import tqdm
from flarelist_generate import estimate_flare_locations, merge_and_process_data
import pandas as pd

# Set to True to fetch new raw_fits, False to use existing CSV
should_fetch_flare_list = False
if __name__ == "__main__":
    tstart = Time("2024-10-01")
    tend = Time("2024-10-02")
    
    if should_fetch_flare_list:
        # step 1
        print("Fetching new flare raw_fits...")
        flare_list = fetch_operational_flare_list(tstart, tend)
        print(f"Found {len(flare_list)} flares in operational list")

        # step 2
        print("Filtering flares and associating CPD files...")
        flare_list_with_files = filter_and_associate_files(flare_list, 'test')
        flare_list_with_files.to_csv('test_flare_data.csv', index=False)
        print(f"Found {len(flare_list_with_files)} flares with CPD files")
    else:
        print("Loading existing flare raw_fits...")
        flare_list_with_files = pd.read_csv('test_flare_data.csv')
        print(f"Loaded {len(flare_list_with_files)} flares from CSV")

    # step 3 calculate and append raw counts

    print("gettin raw counts")

    flare_list_with_files_raw_counts = add_raw_counts_data(flare_list_with_files, save_csv=True)


    # step 4: estimate flare locations and calculate meta pixels and raw counts internally to test
    print("Processing flares with location estimation...")
    flare_list_with_locations = estimate_flare_locations(flare_list_with_files_raw_counts, save_csv=True)

    # step 5: get more coordinate information and tidy
    # final_flarelist_with_locations = merge_and_process_data(flare_list_with_locations)

    print("Processing complete!")




