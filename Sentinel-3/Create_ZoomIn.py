import time

start_time = time.time()

from Sentinel3Class import sentinel3
from Sentinel3_Utilities import get_specific_date_files, create_zoom_images

####################################
##         LOAD CASE INPUTS       ##
####################################

from Inputs import algeria_accident_input, iraq_input, kazakhstan_input, moscow_1_input, moscow_2_input, moscow_both_input, permian_input, varon_algeria_input, varon_turkmenistan_input

input_dict = {
#     'Algeria_Accident': algeria_accident_input,
#     'Permian': permian_input,
#     'Moscow': moscow_both_input,
#     'Moscow_1': moscow_1_input,
    'Moscow_2': moscow_2_input,
#     'Varon_Algeria': varon_algeria_input,
#     'Varon_Turkmenistan': varon_turkmenistan_input
        }

########################################

# Loop over the locations
for key in input_dict.keys():
    
    plume_name = input_dict[key]['plume_name']
    cen_lon = input_dict[key]['cen_lon']
    cen_lat = input_dict[key]['cen_lat']
    all_data_dir = input_dict[key]['all_data_dir']
    specific_main_dates = input_dict[key]['specific_main_dates']    
    main_file_list = get_specific_date_files(specific_main_dates, all_data_dir)
    
    # Loop over main date files of current location
    for main_file in main_file_list[:1]:
        S3_object = sentinel3(plume_name,
                              all_data_dir,
                              main_file,
                              cen_lon,
                              cen_lat)

        # Create ZoomIn Images
        create_zoom_images(plume_name,
                           cen_lat,
                           cen_lon,
                           S3_object.main_date
                          )

########################################################################
# END
########################################################################
    
end_time = time.time()
print()
print('###')
print()
print('Script duration: ', (end_time-start_time), ' seconds')