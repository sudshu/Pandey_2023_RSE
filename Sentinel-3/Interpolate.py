# Interpolate all files to a regular grid
import os
import numpy as np
import numpy.ma as ma

from Sentinel3_Utilities import read_full_slstr_data, get_data_on_grid_of_interest, perform_interpolation, \
                            set_nan_inf_neg_to_zero, get_inland_water_mask

####################################
##              Inputs            ##
####################################

cen_lon = 38.66
cen_lat = 56.61

delta_lonlat = 0.3
pixel_size = 0.001

case_name = 'Moscow'

original_file_dir = '/data/maartenvn/Filtered_data/'+ case_name + '/'

interpolated_dir = '/data/maartenvn/Filtered_data/' + case_name + '_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'

if not os.path.isdir(interpolated_dir):
    os.makedirs(interpolated_dir)

files_to_interpolate = os.listdir(original_file_dir)
files_to_interpolate.sort()
print(delta_lonlat)
####################################
##           INTERPOLATE          ##
####################################

# Define area of interest
grid_of_interest = [cen_lon-delta_lonlat, cen_lon+delta_lonlat, cen_lat-delta_lonlat, cen_lat+delta_lonlat]  

# Desired grid for interpolation later on
lon_reg = np.arange(grid_of_interest[0], grid_of_interest[1]+pixel_size, pixel_size)
lat_reg = np.arange(grid_of_interest[2], grid_of_interest[3]+pixel_size, pixel_size)
lon_reg, lat_reg = np.meshgrid(lon_reg, lat_reg)

print('Interpolating ', len(files_to_interpolate), ' files ... ')

for file in files_to_interpolate:
    lon, lat, s6, s5, s4, s3, s2, s1, landscape, sza, vza, lon_tx, lat_tx = read_full_slstr_data(original_file_dir + file)
    lon, lat, [s6, s5, s4, s3, s2, s1, landscape] = get_data_on_grid_of_interest(lon, lat, [s6, s5, s4, s3, s2, s1, landscape], grid_of_interest)
    water_mask = get_inland_water_mask(landscape)
    lon_tx, lat_tx, [sza, vza] = get_data_on_grid_of_interest(lon_tx, lat_tx, [sza, vza], grid_of_interest)
    sza_average = np.mean(sza)
    vza_average = np.mean(vza)
        
    try:
        [s6, s5, s4, s3, s2, s1, water_mask] = perform_interpolation(lon, lat, lon_reg, lat_reg, [s6, s5, s4, s3, s2, s1, water_mask])
    except:
        continue
    
    [s6, s5, s4, s3, s2, s1, water_mask] = set_nan_inf_neg_to_zero([s6, s5, s4, s3, s2, s1, water_mask])

    # Save lon, lat, s6, s5, s4, s3, s2, s1 and water_mask to interpolated directory
    save_to_path = interpolated_dir + file[:-5]
    if not os.path.isdir(save_to_path):
                os.makedirs(save_to_path)
    np.save(save_to_path + '/lon_interpolated.npy', lon_reg)
    np.save(save_to_path + '/lat_interpolated.npy', lat_reg)
    np.save(save_to_path + '/S6_interpolated.npy', s6)
    np.save(save_to_path + '/S5_interpolated.npy', s5)
    np.save(save_to_path + '/S4_interpolated.npy', s4)
    np.save(save_to_path + '/S3_interpolated.npy', s3)
    np.save(save_to_path + '/S2_interpolated.npy', s2)
    np.save(save_to_path + '/S1_interpolated.npy', s1)
    np.save(save_to_path + '/water_interpolated.npy', water_mask)
    np.save(save_to_path + '/SZA_average.npy', sza_average)
    np.save(save_to_path + '/VZA_average.npy', vza_average)
    
print(' ... Finished interpolation')
