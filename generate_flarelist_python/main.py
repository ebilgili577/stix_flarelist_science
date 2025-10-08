from flarelist_generate import get_flares
from astropy.time import Time
from flarelist_generate import fetch_operational_flare_list, filter_and_associate_files, estimate_flare_locations_and_attenuator, merge_and_process_data

if __name__ == "__main__":
    # manual pipeline
    tstart = Time("2024-10-01")
    tend = Time("2024-10-02")



    # step 1: Fetch the operational flare list
    flare_list = fetch_operational_flare_list(tstart, tend)

    # step 2: filter to counts about 100 and get list of cpd files associated with each
    flare_list_with_files = filter_and_associate_files(flare_list, 'test')

    # step 3: estimate flare locations and get attenuator status
    flare_list_with_locations = estimate_flare_locations_and_attenuator(flare_list_with_files)

    # step 4: get more coordinate information and tidy
    final_flarelist_with_locations = merge_and_process_data(flare_list_with_locations)


    final_flarelist_with_locations[:1]

   