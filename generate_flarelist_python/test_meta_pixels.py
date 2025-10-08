from flarelist_generate import fetch_operational_flare_list, filter_and_associate_files
from astropy.time import Time



from flarelist_generate import estimate_flare_locations_and_attenuator, merge_and_process_data

if __name__ == "__main__":



    tstart = Time("2024-10-01")
    tend = Time("2024-10-02")
    
    flare_list = fetch_operational_flare_list(tstart, tend)

    flare_list_with_files = filter_and_associate_files(flare_list, 'test')

    flare_list_with_files.to_csv('test_flare_data.csv', index=False)

    # flare_list_with_files = pd.read_csv('test_flare_data.csv')

    # step 3: estimate flare locations and calculate meta pixels and raw counts internally to test
    flare_list_with_locations = estimate_flare_locations_and_attenuator(flare_list_with_files, save_csv=True)

