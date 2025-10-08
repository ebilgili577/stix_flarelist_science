from flarelist_generate import fetch_operational_flare_list, filter_and_associate_files
from create_raw_meta_pixels import create_raw_meta_pixels
from astropy.time import Time
from astropy import units as u
from stixpy.product import Product
import numpy as np

from sunpy.time import parse_time


from flarelist_generate import estimate_flare_locations_and_attenuator, merge_and_process_data
import pandas as pd

if __name__ == "__main__":
   # print("Testing create_raw_meta_pixels method directly...")

    # Step 1 & 2: Get flare data (run once, then comment out to avoid re-downloading)
    tstart = Time("2024-10-01")
    tend = Time("2024-10-02")
    
    #print("Fetching operational flare list...")
    #flare_list = fetch_operational_flare_list(tstart, tend)
    #print(f"Found {len(flare_list)} flares")
    
    #print("Filtering flares and associating CPD files...")
    #flare_list_with_files = filter_and_associate_files(flare_list, 'test')
    #print(f"Found {len(flare_list_with_files)} flares with CPD files")
    
    # Save this data so you don't need to re-run steps 1-2
    #flare_list_with_files.to_csv('test_flare_data.csv', index=False)
    #print("Saved flare data to 'test_flare_data.csv' - you can comment out steps 1-2 next time")

    flare_list_with_files = pd.read_csv('test_flare_data.csv')
    print(f"Found {len(flare_list_with_files)} flares with CPD files")
      # step 3: estimate flare locations and get attenuator status
    #flare_list_with_locations = estimate_flare_locations_and_attenuator(flare_list_with_files)



# later
    # step 4: get more coordinate information and tidy
    #final_flarelist_with_locations = merge_and_process_data(flare_list_with_locations)


    results = {"loc_x": [], "loc_y": [], "loc_x_stix": [], "loc_y_stix": [],
               "sidelobes_ratio": [], "flare_id": [], "error": [], "attenuator": []}

    for i, row in flare_list_with_files.iterrows():
        energy_range = [4, 16] * u.keV

        # Define a 20s time range around peak time
        tstart = parse_time(row["peak_UTC"]) - 20 * u.s
        tend = parse_time(row["peak_UTC"]) + 20 * u.s
        time_range = [tstart.strftime("%Y-%m-%dT%H:%M:%S"), tend.strftime("%Y-%m-%dT%H:%M:%S")]
        cpd_file = row["filenames"]
        att = False  # Default value for attenuator

        try:
            cpd_sci = Product(cpd_file)

            # Check for attenuator status by looking for any 'rcr' data points in the time range
            # as the att_in column in the operational flarelist isnt working.
            if np.any(cpd_sci.data[(cpd_sci.data["time"] >= tstart) & (cpd_sci.data["time"] <= tend)]["rcr"]):
                att = True
                energy_range = [4, 25] * u.keV
            print(att)
            

            print('estimating flare location')
            # Estimate flare location

            ## change
            meta_pixels_sci = create_raw_meta_pixels(cpd_sci, 
                                         time_range=time_range, 
                                         energy_range=energy_range, 
                                         flare_location=[0, 0] * u.arcsec, 
                                         no_shadowing=True,
                                         return_raw_counts=True)

        except Exception as e:
            results["loc_x"].append(np.nan)
            results["loc_y"].append(np.nan)
            results["loc_x_stix"].append(np.nan)
            results["loc_y_stix"].append(np.nan)
            results["sidelobes_ratio"].append(np.nan)
            results["error"].append(True)
            results["flare_id"].append(row["flare_id"])
            results["attenuator"].append(att)