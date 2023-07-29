import os, time
start_time = time.time()

from Sentinel3Class import sentinel3
from Sentinel3_Utilities import get_specific_date_files

####################################
##   DEFINE PROCESSING SETTINGS   ##
####################################

##### Reference date selection criteria -> optional parameters when initializing sentinel3() object
min_vis_corr = 0.3
min_delR = -5.
number_of_references_original = 5
reference_day_range = 30

##### Filter settings

destripe = True

filter_delR_mask = False
delR_mask_limit = 0.

filter_auxiliary_mask = True
correlation_minimum = 0.5
auxiliary_bands_to_consider = [1,2,3] # 0=S4, 1=S3, 2=S2, 3=S1

filter_top_xperc = False
top_xperc = 0

filter_median = False
median_filter_size = 0

####################################
##      SET DESIRED RESULTS       ##
####################################

SHOW_STANDARD_TROPOMI = True
SHOW_STANDARD_S3 = True
SHOW_STANDARD_S2 = True

save_results = False

img_quality_S3= 10

destripe_tropomi = True
delta_lonlat_tropomi = 2.0

search_ref_date_s2 = True
multi_refs_S2 = True
area_S2 = 4

VERBOSE = True

####################################
##         LOAD CASE INPUTS       ##
####################################

from Inputs import algeria_accident_input, iraq_input, kazakhstan_input, moscow_1_input, moscow_2_input, moscow_both_input, permian_input, varon_algeria_input, varon_turkmenistan_input, iraq_input, kazakhstan_input

input_dict = {
    'Algeria_Accident': algeria_accident_input,
#     'Iraq': iraq_input,
#     'Kazakhstan': kazakhstan_input,
#     'Permian': permian_input,
#     'Moscow': moscow_both_input,
#     'Moscow_1': moscow_1_input,
#     'Moscow_2': moscow_2_input,
#     'Varon_Algeria': varon_algeria_input,
#     'Varon_Turkmenistan': varon_turkmenistan_input
        }


########################################

# Loop over the locations
for key in input_dict.keys():
    
    if VERBOSE:
        print();print('###################################################');print()
        print(key)
        print();print('###################################################');print()
                
    plume_name = input_dict[key]['plume_name']
    cen_lon = input_dict[key]['cen_lon']
    cen_lat = input_dict[key]['cen_lat']
    delta_lonlat = input_dict[key]['delta_lonlat']
    pixel_size = input_dict[key]['pixel_size']
    all_data_dir = input_dict[key]['all_data_dir']
    specific_main_dates = input_dict[key]['specific_main_dates']
    previous_main_date = None

    if VERBOSE:
        print(all_data_dir)
        
    if specific_main_dates == 'All':
        if VERBOSE:
            print('Processing all dates')
        main_file_list = os.listdir(all_data_dir)
    else:
        if VERBOSE:
            print('Processing specific dates only')
        main_file_list = get_specific_date_files(specific_main_dates, all_data_dir)
        
    # Loop over main date files of current location
    for main_file in main_file_list[:1]:
        
        # Process Sentinel-3 data
        S3_object = sentinel3(plume_name,
                              all_data_dir,
                              main_file,
                              cen_lon,
                              cen_lat,
                              delta_lonlat=delta_lonlat,
                              pixel_size=pixel_size,
                              verbose=VERBOSE,
                              number_of_references_original=number_of_references_original,
                              reference_day_range=reference_day_range
                             )
            
        S3_object.get_correlation_based_reference_data()
        S3_object.perform_multi_band_multi_pass_method()
        S3_object.apply_filters_and_calculate_omega(destripe=destripe,
                                                    filter_delR_mask=filter_delR_mask,
                                                    delR_mask_limit=delR_mask_limit,
                                                    filter_auxiliary_mask=filter_auxiliary_mask,
                                                    correlation_minimum=correlation_minimum,
                                                    auxiliary_bands_to_consider=auxiliary_bands_to_consider,
                                                    filter_top_xperc=filter_top_xperc,
                                                    top_xperc=top_xperc,
                                                    filter_median=filter_median,
                                                    median_filter_size=median_filter_size
                                                   )

        # Process TROPOMI and Sentinel-2, save data and create standard images
        if S3_object.main_date.date() == previous_main_date:
            S3_object.create_standard_images_and_save_data_to_pickle_files(S3_img_quality=img_quality_S3,
                                                                           S3_show=SHOW_STANDARD_S3,
                                                                           retrieve_tropomi=False,
                                                                           retrieve_S2=False,
                                                                           save_results=save_results)

        else:       
            S3_object.create_standard_images_and_save_data_to_pickle_files(S3_img_quality=img_quality_S3,
                                                                           S3_show=SHOW_STANDARD_S3,
                                                                           retrieve_tropomi=True,
                                                                           tropomi_delta_lonlat=delta_lonlat_tropomi,
                                                                           tropomi_show=SHOW_STANDARD_TROPOMI,
                                                                           tropomi_destripe=destripe_tropomi,
                                                                           retrieve_S2=True,
                                                                           S2_area=area_S2,
                                                                           S2_show=SHOW_STANDARD_S2,
                                                                           S2_search_ref_date=search_ref_date_s2,
                                                                           S2_multi_refs=multi_refs_S2,
                                                                           save_results=save_results)
            
            
        previous_main_date = S3_object.main_date.date()

    
########################################################################
# END
########################################################################
    
end_time = time.time()
print()
print('###')
print()
print('Script duration: ', (end_time-start_time), ' seconds')