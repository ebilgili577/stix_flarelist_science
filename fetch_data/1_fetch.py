import sys
import os

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), 'generate_flarelist_python'))

from generate_flarelist_python.flarelist_generate import *

tstart = '2021-01-01'
tend = '2021-12-31'

tstart = Time(tstart)
tend = Time(tend)

# path where fits files are stored locally
raw_fits_path = 'raw_fits'

# step 1: Fetch the operational flare list
flare_list = fetch_operational_flare_list(tstart, tend)

# step 2: filter to counts about 100 and get list of cpd files associated with each
flare_list_with_files = filter_and_associate_files(flare_list, raw_fits_path, save_csv=True)

# step 3: add raw counts and attenuator status to the list
flare_list_with_files_raw_counts = add_raw_counts_data(flare_list_with_files, save_csv=True)

# step 4: estimate flare locations
flare_list_with_locations = estimate_flare_locations(flare_list_with_files_raw_counts, save_csv=True)

# step 5: get more coordinate information and tidy
final_flarelist_with_locations = merge_and_process_data(flare_list_with_locations, save_csv=True)