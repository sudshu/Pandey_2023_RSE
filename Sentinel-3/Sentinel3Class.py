import sys, glob, math
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td

from scipy import ndimage, fftpack
from skimage.metrics import structural_similarity

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy import geodesic
from shapely.geometry import LineString
import cv2 as cv
import pickle

sys.path.append("./Code_Tropomi")
sys.path.append("./Code_Tropomi/helper")
sys.path.append("./Code_Sentinel_2")

from Sentinel3_Utilities import read_full_interpolated_data, normalize, to_integer, image_spoof, get_data_on_grid_of_interest

# Import Tropomi data processing
from TROPOMI_orbit_plotter import makeOrbitWisePlotsCartopy
from region_rectangle_defination import makeCordbox

# Import Sentinel 2 data processing
import MS_XCH4_retrieval as ms


class sentinel3:
    """
    
    Object used to analyze the band S6 (SWIR-2) data from the Sentinel 3 satellites (A and B) in order to identify and quantify methane super emissions
    
    """
    
    def __init__(self,
                 plume_name: str,
                 file_dir: str,
                 main_date_file: str,
                 lon: float,
                 lat: float,
                 delta_lonlat: float=0.15,
                 pixel_size: float=0.001,
                 verbose: bool=False,
                 min_vis_corr: float=-1.,
                 min_delR: float = -5.0,
                 number_of_references_original: int=5,
                 reference_day_range: int = 20
                ):
        
        """
        Intitialize the Sentinel 3 object
        
        Inputs:
        - plume_name: str                        name of location or plume
        - file_dir: str                          directory where S3 files are stores
        - file: str                              file to be analyzed
        - lon: float                             central longitude of location to be analyzed
        - lat: float                             central latitude of location to be analyzed
        - delta_lonlat: float                    half width/height of final image in degrees
        - minimum_delR_correlation: float        minimum delR correlation between main and reference
        - number_of_references_original: int     number of reference dates that the median will be taken of for the final reference data
        - minimum_delR: float                    lowest allowed delR MBMP value, smaller values are assumed to not be realistic
        - reference_day_range: int               number of days around the main date where for reference dates is searched
        - verbose: bool                          print information while processing

        Defines in addition:
        - self.main_date
        - self.ref_file_list
        - self.grid_of_interest
        - self.lon_grid
        - self.lat_grid

        """
        
        self.plume_name = plume_name
        self.file_dir = file_dir
        self.main_date_file = main_date_file
        self.lon = lon
        self.lat = lat
        self.delta_lonlat = delta_lonlat
        self.pixel_size = pixel_size
        self.verbose = verbose
        self.min_vis_corr = min_vis_corr
        self.min_delR = min_delR
        self.number_of_references_original = number_of_references_original
        self.reference_day_range = reference_day_range
        self.valid = True
                
        # Get main date datetime object from file name
        self.main_date = dt.strptime(self.main_date_file.split('/')[-1][16:31], '%Y%m%dT%H%M%S')
        self.main_satellite = str(self.main_date_file.split('/')[-1][0:3])
        
        self.swir2_main, self.swir1_main, s4_main, s3_main, s2_main, s1_main, self.water_main_mask, self.sza_main, self.vza_main = read_full_interpolated_data(self.main_date_file)
        self.auxiliary_bands_main = [s4_main, s3_main, s2_main, s1_main]    
            
        if self.verbose:
            print('\n \n ### \n \n', self.main_date, '\n')
            
        # Get all available files
        all_files_list = glob.glob(self.file_dir + '*')
        
        # Get main date files (should have 1 or 2 items)
        all_main_date_files = glob.glob(self.file_dir + '*__' + self.main_date.strftime("%Y%m%d") + 'T*')
        
        # Get reference date files, excluding the main date files
        main_date_file_indices = []
        for main_date_file in all_main_date_files:
            main_date_file_index = list(all_files_list).index(main_date_file)
            all_files_list = np.concatenate((all_files_list[:main_date_file_index], all_files_list[main_date_file_index+1:]))
        self.ref_file_list = all_files_list      

        # Define lonlat grid
        self.grid_of_interest = [self.lon-self.delta_lonlat, self.lon+self.delta_lonlat, self.lat-self.delta_lonlat, self.lat+self.delta_lonlat]  
#         self.grid_of_interest = [self.lon-0.225, self.lon+0.225, self.lat-0.175, self.lat+0.175]   # Used for Moscow two-plume-in-one-image case
        lon_reg = np.arange(self.grid_of_interest[0], self.grid_of_interest[1]+self.pixel_size, self.pixel_size)
        lat_reg = np.arange(self.grid_of_interest[2], self.grid_of_interest[3]+self.pixel_size, self.pixel_size)
        self.lon_grid, self.lat_grid = np.meshgrid(lon_reg, lat_reg)

    def get_correlation_based_reference_data(self):
        
        """
        Searches for and returns valid reference dates within the main file directory based on the correlation of the visual S1 band and the desired number of references
        
        Returns:
        - ref_file_to_use: list of str
        
        """
        corr_list = []
        ref_files_considered = []
        
        # Obtain main date visual observations
        s1_main = np.load(self.main_date_file + '/S1_interpolated.npy')
        
        dates_within_range = 0
        for ref_file in self.ref_file_list:
            # Get reference datetime object from file name            
            ref_date = dt.strptime(ref_file.split('/')[-1][16:31], '%Y%m%dT%H%M%S')
            
            # Only consider the date if within the desired timeframe
            date_difference = abs((self.main_date - ref_date).days)
            if date_difference <= self.reference_day_range:
                dates_within_range += 1
                
                # Read reference date data and calculate visual band correlation
                s1_ref = np.load(ref_file + '/S1_interpolated.npy')
                correlation = np.corrcoef(s1_main.flatten(), s1_ref.flatten())[0,1]
                if not math.isnan(correlation):
                    ref_files_considered.append(ref_file)
                    corr_list.append(correlation)
                    
        # Pick the highest correlation dates as references to use
        number_of_references = self.number_of_references_original            
        if len(corr_list) >= number_of_references:
            indices_to_use = np.argpartition(corr_list, -number_of_references)[-number_of_references:] 
            ref_files_to_use = [ref_files_considered[indices_to_use[i]] for i in range(number_of_references)]
            best_correlations = [corr_list[indices_to_use[i]] for i in range(number_of_references)]
            if self.verbose:
                print('Band S1 correlations of best ', number_of_references, ' ref dates: ', best_correlations) 
        elif len(corr_list) > 0:
            number_of_references = len(corr_list)
            indices_to_use = np.argpartition(corr_list, -number_of_references)[-number_of_references:] 
            ref_files_to_use = [ref_files_considered[indices_to_use[i]] for i in range(number_of_references)]
            best_correlations = [corr_list[indices_to_use[i]] for i in range(number_of_references)]
            if self.verbose:
                print('Band S1 correlations of best ', number_of_references, ' ref dates: ')
                print(best_correlations) 
        else:
            if self.verbose:
                print('No reference dates within ', self.reference_day_range, ' days available / good enough ->  discarded')
                print('Number of dates within valid range, but NaN correlation: ', dates_within_range)
            ref_files_to_use = None
            best_correlations = None
            self.valid = False
                
        self.ref_files_to_use = ref_files_to_use

    
    def perform_multi_band_multi_pass_method(self):
        
        """
        
        Process the data according to the Multi-Band Multi-Pass method (Varon et al., 2021). As reference data, the median of multiple dates is taken (as in ref_files_to_use)
        It also retrieves the corresponding wind data from ERA5.
                
        """
                      
        if self.valid:
            self.number_of_references = len(self.ref_files_to_use)

            if self.verbose:
                print('Reference dates:')
            
            if self.number_of_references > 1:
                self.ref_satellite = str(self.number_of_references) + ' refs'
                self.ref_date = str(self.number_of_references) + ' refs'
                self.sza_ref = str(self.number_of_references) + ' refs'
                self.vza_ref = str(self.number_of_references) + ' refs'
                
                all_s6_ref = []; all_s5_ref = []; all_s4_ref = []; all_s3_ref = []; all_s2_ref = []
                all_s1_ref = []; all_water_mask_ref = []; all_delR_mbmp = []; all_delR_ref = []

                # Perform delR MBMP on each of the reference dates
                for ref_file in self.ref_files_to_use:
                    current_reference_date = dt.strptime(ref_file.split('/')[-1][16:31], '%Y%m%dT%H%M%S')

                    if self.verbose:
                        print('       ', current_reference_date)
                    
                    swir2_ref, swir1_ref, s4_ref, s3_ref, s2_ref, s1_ref, water_ref_mask, sza_ref, vza_ref = read_full_interpolated_data(ref_file)                   
                    delR_mbmp_original, delR_main_original, delR_ref_original = self.calculate_delR(swir1_ref, swir2_ref)
    
                    all_delR_mbmp.append(delR_mbmp_original); all_delR_ref.append(delR_ref_original)
                    all_s6_ref.append(swir2_ref); all_s5_ref.append(swir1_ref); all_s4_ref.append(s4_ref) 
                    all_s3_ref.append(s3_ref); all_s2_ref.append(s2_ref); all_s1_ref.append(s1_ref) 
                    all_water_mask_ref.append(water_ref_mask)

                # Take median to form one reference observation
                self.swir2_ref = np.median(np.array(all_s6_ref), axis=0)
                self.swir1_ref = np.median(np.array(all_s5_ref), axis=0)   
                s4_ref = np.median(np.array(all_s4_ref), axis=0)
                s3_ref = np.median(np.array(all_s3_ref), axis=0)   
                s2_ref = np.median(np.array(all_s2_ref), axis=0)
                s1_ref = np.median(np.array(all_s1_ref), axis=0)
                self.water_ref_mask = np.median(np.array(all_water_mask_ref), axis=0)
                self.auxiliary_bands_ref = [s4_ref, s3_ref, s2_ref, s1_ref]

                self.delR_ref_original = np.median(np.array(all_delR_ref), axis=0)
                self.delR_main_original = delR_main_original
                
                # Recalculate delR MBMP based on median of delR MBSP for multiple reference dates
                self.delR_mbmp_original = self.delR_main_original - self.delR_ref_original                 

            else: 
                ref_file = self.ref_files_to_use[0]
                self.ref_satellite = str(ref_file.split('/')[-1][0:3])
                self.ref_date = dt.strptime(ref_file.split('/')[-1][16:31], '%Y%m%dT%H%M%S')

                self.swir2_ref, self.swir1_ref, s4_ref, s3_ref, s2_ref, s1_ref, self.water_ref_mask, self.sza_ref, self.vza_ref = read_full_interpolated_data(ref_file)
                self.auxiliary_bands_ref = [s4_ref, s3_ref, s2_ref, s1_ref]

                self.delR_mbmp_original, self.delR_main_original, self.delR_ref_original = self.calculate_delR(self.s5_ref, self.s6_ref)
                
            # Calculate the Pearson product-moment correlation coefficients between main date and new reference observations
            self.delR_mbsp_corr = np.corrcoef(self.delR_main_original.flatten(), self.delR_ref_original.flatten())[0,1]
            self.aux_corr = np.corrcoef(self.auxiliary_bands_main[-1].flatten(), self.auxiliary_bands_ref[-1].flatten())[0,1]
            
            if self.verbose:
                print('New Band S1 correlation: ', self.aux_corr)
            
            if self.aux_corr < self.min_vis_corr: 
                if self.verbose:
                    print('New Band S1 correlation below ', self.min_vis_corr, ' -> discarded')
                self.valid = False

            #Get wind data through the Sentinel-2 code
            s2_obj_main = ms.ch4ret(self.main_date, self.main_date,'MBMP', lat=self.lat, lng=self.lon, area= 3)
            s2_obj_main_minus_1hr = ms.ch4ret(self.main_date-td(hours=1), self.main_date-td(hours=1), 'MBMP', lat=self.lat, lng=self.lon, area=3)
            s2_obj_main_minus_2hr = ms.ch4ret(self.main_date-td(hours=2), self.main_date-td(hours=2), 'MBMP', lat=self.lat, lng=self.lon, area=3)

            s2_obj_main.getWindSpeed()
            s2_obj_main_minus_1hr.getWindSpeed()
            s2_obj_main_minus_2hr.getWindSpeed()

            u10_main, v10_main, self.wind_main, self.wind_dir_main = np.array([s2_obj_main.u10]), np.array([s2_obj_main.v10]), s2_obj_main.wind_speed, s2_obj_main.wind_dir   
            u10_main_minus_1hr, v10_main_minus_1hr, self.wind_main_minus_1hr, self.wind_dir_main_minus_1hr = np.array([s2_obj_main_minus_1hr.u10]), np.array([s2_obj_main_minus_1hr.v10]), s2_obj_main_minus_1hr.wind_speed, s2_obj_main_minus_1hr.wind_dir
            u10_main_minus_2hr, v10_main_minus_2hr, self.wind_main_minus_2hr, self.wind_dir_main_minus_2hr = np.array([s2_obj_main_minus_2hr.u10]), np.array([s2_obj_main_minus_2hr.v10]), s2_obj_main_minus_2hr.wind_speed, s2_obj_main_minus_2hr.wind_dir

            self.u10_main_norm, self.v10_main_norm = u10_main/self.wind_main, v10_main/self.wind_main
            self.u10_main_minus_1hr_norm, self.v10_main_minus_1hr_norm = u10_main_minus_1hr/self.wind_main_minus_1hr, v10_main_minus_1hr/self.wind_main_minus_1hr
            self.u10_main_minus_2hr_norm, self.v10_main_minus_2hr_norm = u10_main_minus_2hr/self.wind_main_minus_2hr, v10_main_minus_2hr/self.wind_main_minus_2hr
            
            if type(self.ref_date) != str:
                s2_obj_ref = ms.ch4ret(self.ref_date, self.ref_date, 'MBMP', lat=self.lat, lng=self.lon, area=3)
                s2_obj_ref.getWindSpeed()
                u10_ref, v10_ref, self.wind_ref, self.wind_dir_ref = np.array([s2_obj_ref.u10]), np.array([s2_obj_ref.v10]), s2_obj_ref.wind_speed, s2_obj_ref.wind_dir
                self.u10_ref_norm, self.v10_ref_norm = u10_ref/self.wind_ref, v10_ref/self.wind_ref
            else: 
                self.wind_ref = str(self.number_of_references) + ' refs'
                self.wind_dir_ref = str(self.number_of_references) + ' refs'
                self.u10_ref_norm = 0.
                self.v10_ref_norm = 0.
               
        else:
            self.valid = False


    def apply_filters_and_calculate_omega(self,
                                          destripe=False,
                                          filter_delR_mask=True,
                                          delR_mask_limit=0.,
                                          filter_auxiliary_mask = True,
                                          correlation_minimum = 0.5,
                                          auxiliary_bands_to_consider = [1,2,3], # 0=S4, 1=S3, 2=S2, 3=S1
                                          filter_top_xperc = False,
                                          top_xperc = 0,
                                          filter_median = False,
                                          median_filter_size = 0
                                         ):

        """
        
        Generates plume filters on the original delR MBMP array. The first filter is based on the structural similarity between the delR MBSP images of main and reference dates, to identify regions of interest, where potential plumes exist. The second filter is based on the structural similarity between the auxiliary bands (S1, S2, S3 and/or S4) of the main and reference dates, to identify regions where characteristics are due to clouds or surface features, where differences in signal are likely NOT due to methane. Lastly a top x percentile filter and median filter can be applied, by default these are turned off.
        
        """
        self.filter_delR_mask = filter_delR_mask
        self.delR_mask_limit = delR_mask_limit
        self.filter_auxiliary_mask = filter_auxiliary_mask
        self.correlation_minimum = correlation_minimum
        self.auxiliary_bands_to_consider = auxiliary_bands_to_consider
        self.filter_top_xperc = filter_top_xperc
        self.top_xperc = top_xperc
        self.filter_median = filter_median
        self.median_filter_size = median_filter_size
        
        if self.verbose:
            print('SZA main =  ', self.sza_main)
            print('VZA main =  ', self.vza_main)
        
        if self.valid:

            # Apply median filter to the inland water mask, because there are a lot of tiny dots on the image due to this mask
            self.water_main_mask = ndimage.median_filter(self.water_main_mask, size=2)
            self.water_ref_mask = ndimage.median_filter(self.water_ref_mask, size=2)
            
            # Mask delR MBMP with inland water mask
            self.delR_mbmp_no_water = self.delR_mbmp_original.copy()
            self.delR_mbmp_no_water[np.where( self.water_main_mask < 1)] = 0
            self.delR_mbmp_no_water[np.where( self.water_ref_mask < 1)] = 0
                    
            # Apply destriping procedure
            if destripe:
                # FFT of stripes enhanced image
                delR_mbmp_no_water_blurred = ndimage.gaussian_filter(self.delR_mbmp_no_water, 1.)
                stripes = self.delR_mbmp_no_water - delR_mbmp_no_water_blurred
                stripes_fft = fftpack.fft2(stripes)
                stripes_fft_centered = np.fft.fftshift(stripes_fft)

                # Obtain angle from FFT
                stripes_fft_max_index = np.argwhere(np.abs(stripes_fft_centered) == np.nanmax(np.abs(stripes_fft_centered)))
                x_range = (stripes_fft_max_index[1][1] - stripes_fft_max_index[0][1])
                y_range = (stripes_fft_max_index[1][0] - stripes_fft_max_index[0][0])
                perp_angle = np.rad2deg(np.arctan2(y_range, x_range))
                angle = perp_angle-90
                
                if self.verbose:
                    print('Stripe angle: ', angle, ' degrees')
                
                # Subtract median per row / stripe (after rotation)
                delR_rotated = ndimage.rotate(self.delR_mbmp_no_water, angle, reshape=False)
                median_per_row = np.zeros(delR_rotated.shape)
                delR_minus_median = delR_rotated.copy()
                for i in range(len(delR_rotated)):
                    median_per_row[i] = np.nanmedian(delR_rotated[i])
                    delR_minus_median[i] -= median_per_row[i]
                self.delR_mbmp_destriped = ndimage.rotate(delR_minus_median, -angle, reshape=False)
                
                # Fill values that were cut in the destriping processs to their original values
                for i in range(len(self.delR_mbmp_destriped)):
                    for j in range(len(self.delR_mbmp_destriped[0])):
                        if self.delR_mbmp_destriped[i][j] == 0.:
                            self.delR_mbmp_destriped[i][j] = self.delR_mbmp_no_water[i][j]
            else:
                self.delR_mbmp_destriped = self.delR_mbmp_no_water
                            
            # Remove positives
            self.delR_mbmp_limited = self.delR_mbmp_destriped.copy()
            self.delR_mbmp_limited[np.where(self.delR_mbmp_destriped > 0)] = 0

            if self.verbose:
                print('Minimum delR MBMP: ', self.delR_mbmp_limited.min())

            if self.delR_mbmp_limited.min() < self.min_delR:
                if self.verbose:
                    print('Unrealistically low minimum delR MBMP ', self.delR_mbmp_limited.min(), '-> discarded')
                self.valid = False
                
            else:
                # Get regions of interest based on delR mask (delR MBSP of main and ref date)
                if self.filter_delR_mask:
                    self.delR_mask, delR_mbmp_interest = self.apply_delR_mask(self.delR_mbmp_limited, self.delR_main_original, self.delR_ref_original, self.delR_mask_limit)
                else:
                    self.delR_mask = 255*np.zeros(self.delR_mbmp_limited.shape)
                    delR_mbmp_interest = self.delR_mbmp_limited.copy()
                self.delR_mask = normalize(self.delR_mask)

                # Mask based on auxiliary bands (noisy features)
                if self.filter_auxiliary_mask:
                    self.auxiliary_mask, delR_mbmp_masked, self.bands_in_aux_mask, correlations_of_mask, scores_of_mask = self.apply_auxiliary_mask(delR_mbmp_interest, self.auxiliary_bands_main, self.auxiliary_bands_ref, self.auxiliary_bands_to_consider, self.correlation_minimum)
                    plt.show()
                else:
                    self.auxiliary_mask = np.zeros(delR_mbmp_interest.shape)
                    delR_mbmp_masked = delR_mbmp_interest.copy()
                self.auxiliary_mask = normalize(self.auxiliary_mask)
                
                # Filter top 100-x% - only remain most strong x% signal
                if self.filter_top_xperc:
                    delR_mbmp_xperc_filtered = delR_mbmp_masked.copy()
                    percentile = np.percentile(delR_mbmp_xperc_filtered.flatten(), self.top_xperc)
                    delR_mbmp_xperc_filtered[delR_mbmp_xperc_filtered > percentile]=0
                else:
                    delR_mbmp_xperc_filtered = delR_mbmp_masked.copy()

                # Median fiLter
                if self.filter_median:
                    delR_mbmp_median_filtered = ndimage.median_filter(delR_mbmp_xperc_filtered, size=self.median_filter_size)
                else:
                    delR_mbmp_median_filtered = delR_mbmp_xperc_filtered.copy()
                
                # Combine all masks into one to save
                self.mask_to_save = self.auxiliary_mask.copy()
                self.mask_to_save[np.where(self.water_main_mask < 1)] = 255
                self.mask_to_save[np.where(self.water_ref_mask < 1)] = 255
                    
                # Apply filter to original delR MBMP
                self.delR_mbmp_final = delR_mbmp_median_filtered.copy()
                self.delR_mbmp_final[np.where(delR_mbmp_median_filtered == 0)] = 0
                
                self.delR_mbmp_final_plume_only = delR_mbmp_median_filtered.copy()
                self.delR_mbmp_final_plume_only[np.where(delR_mbmp_median_filtered == 0)] = 'nan'
                
                # Convert delR (radiance difference) to Omega (methane concentration)
                self.omega = self.fullMBMP2Omega_(self.delR_mbmp_destriped, 'S3', self.sza_main, vza=self.vza_main)
                self.omega_filtered = self.fullMBMP2Omega_(self.delR_mbmp_final_plume_only, 'S3', self.sza_main, vza=self.vza_main)

    
    def calculate_delR(self, swir1_ref, swir2_ref, method='A'):
        
        """
        Calculates the delR arrays using the Multi-Band Multi-Pass method, which uses the delR Multi-Band Single-Pass arrays as an intermediate step.
        
        Input:
        - swir1 band data reference date
        - swir2 band data reference date
        
        Output: 
        - delR Multi-Band Multi-Pass
        - delR Multi-Band Single-Pass main date
        - delR Multi-Band Single-Pass reference  date
        
        """
        
        np.seterr(divide='ignore', invalid='ignore')
        if method=='A':
            delR, delR_main, delR_ref = self.mbmpA(swir2_ref, swir1_ref)
        elif method=='B':
            delR, delR_main, delR_ref = self.mbmpB(swir2_ref, swir1_ref)
        else: print('Invalid method')
            
        infmask_main = np.isinf(delR_main)
        nanmask_main = np.isnan(delR_main)
        delR_main[nanmask_main] = 0
        delR_main[infmask_main] = 0
        
        infmask_ref = np.isinf(delR_ref)
        nanmask_ref = np.isnan(delR_ref)
        delR_ref[nanmask_ref] = 0
        delR_ref[infmask_ref] = 0

        infmask = np.isinf(delR)
        nanmask = np.isnan(delR)
        delR[nanmask] = 0
        delR[infmask] = 0

        return delR, delR_main, delR_ref
    
    def mbmpA(self, swir2_ref, swir1_ref):
        'Initial MBMP method, as proposed by Varon'
        swir2_main = self.swir2_main / np.nanmedian(self.swir2_main)
        swir1_main = self.swir1_main / np.nanmedian(self.swir1_main)
        swir2_ref = swir2_ref / np.nanmedian(swir2_ref)
        swir1_ref = swir1_ref / np.nanmedian(swir1_ref)
        delR_main = np.divide((swir2_main-swir1_main), swir1_main)
        delR_ref = np.divide((swir2_ref-swir1_ref), swir1_ref)
        delR = delR_main - delR_ref
        return delR, delR_main, delR_ref
        
    def mbmpB(self, swir2_ref, swir1_ref): 
        'MBMP method adjusted to match the Radiative transfer model, as used for the delR to Omega transition'
        swir2_main = self.swir2_main / np.nanmedian(self.swir2_main)
        swir1_main = self.swir1_main / np.nanmedian(self.swir1_main)
        swir2_ref = swir2_ref / np.nanmedian(swir2_ref)
        swir1_ref = swir1_ref / np.nanmedian(swir1_ref)
        delR1=np.divide((swir2_main-swir2_ref), swir2_ref)   
        delR2=np.divide((swir2_ref-swir1_ref), swir1_ref)
        delR=delR1-delR2
        return delR, delR1, delR2

    def convert_delR_to_omega(self, delr, satellite, sza,  method , vza= 0):
        
        """
        Transforms delR (radiance) to omgea (mol/m^2) according to Varon-m method (Varon et al., 2021)
        
        Inputs:
        - delr: single number or an 1d or 2d array.    
        - Satellite: S2 or S3
        - sza: solar zenith angle 
        - VZA: viewing zenith angle
        - method: MBMP or MBSP or SBMP

        output
        - omega: methane enhancement in mol/m^2
        
        """
        
        mdata= pickle.load(open('amf_mdata_poly_10_delr_to_omega.pkl', 'rb'))
        tamf = (1/np.cos(np.radians(sza)) + 1/np.cos(np.radians(sza)))
        amf= list(mdata[satellite].keys())[np.argmin(abs(np.array(list(mdata[satellite].keys())) - tamf))]
        if self.verbose:
            print ('AMF=', amf)
        pol = np.poly1d(mdata[satellite][amf][method] )
        omega= pol(np.log(delr+1))

        return omega
    
    
    def fullMBMP2Omega_(self, delr, satellite, sza, vza= 0, flag= True ):
        shape_omega=    delr.shape 
        delr= delr.flatten()
        mdata= pickle.load(open('/home/maartenvn/Code/test_srf_210104_delr_to_omega.pkl', 'rb'))
        
        tamf = (1/np.cos(np.radians(sza)) + 1/np.cos(np.radians(sza)))

#         tamf = ms.giveamf(sza,vza)
     #   print ("AMF", tamf, "%2.1f"%tamf)
        ind_0 = np.where(mdata["omegas"]== 0 )
        dat = mdata[satellite]["%2.1f"%tamf]
        omega= np.zeros_like(delr)
        aa= dat["SWIR2"]/ dat["SWIR2"][0]
        bb= dat["SWIR1"]/ dat["SWIR1"][0]   
        if flag:
            mbmp=  aa/bb-1
        else:
            mbmp =  aa- bb
        for ii, drr in enumerate(delr):
            ind = np.argmin(abs(mbmp -  drr))
            omega[ii]=mdata["omegas"][ind]
        omega= omega.reshape(shape_omega)
        
        return omega

    
    def apply_structural_similarity(self, before, after):
        """
        Function to compute the mean structural similarity index between two images
        and to compute the full structural similarity image (diff image).
        A filter is then applied to the diff image to retain only the large differences
        between the images obtaining a mask of the differences.
        input: two images with a single channel (grey scale images)
        output: structural similarity index, structural similarity image, and mask image
        """ 
        # Compute SSIM between two images
        (score, diff) = structural_similarity(before, after, full=True)
        # The diff image contains the actual image differences between the two images
        # and is represented as a floating point data type in the range [0,1]
        # Convert it to 8-bit unsigned integers in the range [0,255]
        diff = to_integer(diff)
        # Threshold the difference image, followed by finding contours to
        # obtain the regions of the two input images that differ
        mask = cv.threshold(diff, 0, 255, cv.THRESH_BINARY_INV | cv.THRESH_OTSU)[1]
        return score, diff, mask
    
    
    def apply_delR_mask(self, unfiltered_mbmp, unfiltered_main, unfiltered_ref, delR_mask_limit):
    
        unfiltered_main_limited = unfiltered_main.copy()
        unfiltered_ref_limited = unfiltered_ref.copy()
        unfiltered_main_limited[np.where(unfiltered_main > delR_mask_limit)] = delR_mask_limit
        unfiltered_ref_limited[np.where(unfiltered_ref > delR_mask_limit)] = delR_mask_limit

        # Apply structural similarity between the original delR main and delR reference to identify noise sources within these bands
        (score, diff, delR_mask) = self.apply_structural_similarity(to_integer(normalize(unfiltered_ref_limited)),
                                                               to_integer(normalize(unfiltered_main_limited)))

        filtered_mbmp = unfiltered_mbmp.copy()
        filtered_mbmp[np.where(delR_mask == 0)] = 0

        return delR_mask, filtered_mbmp
    

    def apply_auxiliary_mask(self, original, bands_main, bands_ref, bands_to_consider, correlation_min):
    
        final_mask = np.zeros(original.shape)
        bands_in_aux_mask = []
        correlations_of_mask = []
        scores_of_mask = []

        for i in bands_to_consider:
            band_main = bands_main[i]
            band_ref = bands_ref[i]
            (score, diff, mask) = self.apply_structural_similarity(to_integer(normalize(band_ref)), to_integer(normalize(band_main)))

            # Only consider the mask if the correlation is high enough
            if score > correlation_min:
                if self.verbose:
                    print('Mask S' + str(4-i) + ' included, correlation=' + str(score) )
                final_mask[np.where(mask == 255)] = 255
                bands_in_aux_mask.append('S' + str(4-i))
                correlations_of_mask.append(np.corrcoef(band_main.flatten(), band_ref.flatten())[0,1])
                scores_of_mask.append(score)

            else:
                if self.verbose:
                    print('Mask S' + str(4-i) + ' NOT included, correlation=' + str(score) )

        delR_mbmp_masked= original.copy()
        delR_mbmp_masked[np.where(final_mask == 255)] = 0

        return final_mask, delR_mbmp_masked, bands_in_aux_mask, correlations_of_mask, scores_of_mask
    
    
    def create_standard_image_and_save_data(self,
                                            img_quality = 10,
                                            show = False,
                                            save = False
                                           ):
        
        if self.valid:
            
            # Filter unrealistically large omega values, as they will influence the colorscale limits
            self.omega[np.where(self.omega > 100000)] = 0

            # Define colorscale limits
            standard_deviations = 3
            vmin_delR = np.mean(self.delR_mbmp_destriped) - standard_deviations * np.std(self.delR_mbmp_destriped)
            vmax_delR = np.mean(self.delR_mbmp_destriped) + standard_deviations * np.std(self.delR_mbmp_destriped)        
            # Get array sizes
            entries_horizontal = self.lon_grid.shape[0]
            entries_vertical = self.lon_grid.shape[1]
            
            # Get grid area in km
            myGeod = geodesic.Geodesic()
            shapelyObject_horizontal = LineString([(self.grid_of_interest[0], self.grid_of_interest[1]), (self.grid_of_interest[2], self.grid_of_interest[2])])
            shapelyObject_vertical = LineString([(self.grid_of_interest[0], self.grid_of_interest[0]), (self.grid_of_interest[2], self.grid_of_interest[3])])
            horizontal_distance = int(myGeod.geometry_length(shapelyObject_horizontal)/1000)
            vertical_distance = int(myGeod.geometry_length(shapelyObject_vertical)/1000)

            # Define coordinates of wind arrows
            lon_wind_1 = np.array([self.lon_grid[int(entries_horizontal / 10)][int(entries_vertical/10)]])
            lat_wind_1 = np.array([self.lat_grid[int(-entries_horizontal / 10)][int(-entries_vertical/10)]])
            lon_wind_2 = np.array([self.lon_grid[int(entries_horizontal / 10)][int(3*entries_vertical/10)]])
            lat_wind_2 = np.array([self.lat_grid[int(-entries_horizontal / 10)][int(-3*entries_vertical/10)]])
            lon_wind_3 = np.array([self.lon_grid[int(entries_horizontal / 10)][int(5*entries_vertical/10)]])
            lat_wind_3 = np.array([self.lat_grid[int(-entries_horizontal / 10)][int(-5*entries_vertical/10)]])

            # Define units for lonlat
            if self.lon >= 0.:
                lon_unit = '$\degree$E'
            else:
                lon_unit = '$\degree$W'
            if self.lat >= 0.:
                lat_unit = '$\degree$N'
            else:
                lat_unit = '$\degree$S'   

            # Load Google Earth background image
            cimgt.QuadtreeTiles.get_image = image_spoof # reformat web request for street map spoofing
            img = cimgt.QuadtreeTiles() # spoofed, downloaded street map

            # Initialize figure
            fig = plt.figure(figsize=(14, 7))
            
            # Figure title
            if type(self.ref_date) != str:
                fig.suptitle(
                    self.plume_name + '  |  ' + str(np.around(self.lon, decimals=2))  + lon_unit +  '  ' + str(np.around(self.lat, decimals=2)) + lat_unit + '  |  ' + 'Area=' + str(horizontal_distance) + 'x' + str(vertical_distance) + '$km^2$ \n' +
                    'Main: ' + self.main_satellite + ' - ' + self.main_date.strftime('%Y-%m-%d %H:%M') +'  |  ' + 'Ref. : ' + self.ref_satellite + ' - ' + self.ref_date.strftime('%Y-%m-%d_%H%M'),
                    fontsize=10                )
            else:
                fig.suptitle(
                    self.plume_name + '  |  ' + str(np.around(self.lon, decimals=2))  + lon_unit +  '  ' + str(np.around(self.lat, decimals=2)) + lat_unit + '  |  ' + 'Area=' + str(horizontal_distance) + 'x' + str(vertical_distance) + '$km^2$ \n' +
                    'Main: ' + self.main_satellite + ' - ' + self.main_date.strftime('%Y-%m-%d %H:%M') + '  |  ' + 'Ref. : ' + self.ref_satellite,
                    fontsize=10
                )
                        
            ####  Original (limited) delR MBMP with wind
            ax = fig.add_subplot(1, 2, 1, projection=img.crs)
            ax.set_title(
                '$\Delta$R MBMP \n' +
                'Wind:   ' + str(np.around(self.wind_main_minus_2hr, decimals=1)) + '    |    ' + str(np.around(self.wind_main_minus_1hr, decimals=1)) + '    |    ' + str(np.around(self.wind_main, decimals=1)) + '   m/s', loc='left',
                fontsize=10)
            ax.set_extent([self.lon-0.15, self.lon+0.15, self.lat-0.15, self.lat+0.15])
#             ax.set_extent([self.lon-0.225, self.lon+0.225, self.lat-0.175, self.lat+0.175])
            ax.coastlines()
            gl = ax.gridlines(draw_labels=False)
            gl.xlabels_top = False
            gl.ylabels_right = False
            pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.delR_mbmp_destriped, cmap=plt.cm.coolwarm, vmin=vmin_delR, vmax=-vmin_delR, shading='auto', transform=ccrs.PlateCarree())
            plt.scatter(self.lon, self.lat, s=100, c='black', marker='x', transform=ccrs.PlateCarree())            
            plt.quiver(lon_wind_1, lat_wind_1, self.u10_main_minus_2hr_norm, self.v10_main_minus_2hr_norm, scale = 10, color='black', transform=ccrs.PlateCarree())
            plt.quiver(lon_wind_2, lat_wind_2, self.u10_main_minus_1hr_norm, self.v10_main_minus_1hr_norm, scale = 10, color='black', transform=ccrs.PlateCarree())
            plt.quiver(lon_wind_3, lat_wind_3, self.u10_main_norm, self.v10_main_norm, scale = 10, color='black', transform=ccrs.PlateCarree())

            ####  Filtered plume figure with background Google Earth image
            ax = fig.add_subplot(2, 4, 3, projection=img.crs)
            ax.set_title('$\Delta$R MBMP   |   $R_{MBSP}$=' + str(np.around(float(self.delR_mbsp_corr), decimals=3)), fontsize=10)
            ax.set_extent([self.lon-0.15, self.lon+0.15, self.lat-0.15, self.lat+0.15], crs=ccrs.PlateCarree())
#             ax.set_extent([self.lon-0.225, self.lon+0.225, self.lat-0.175, self.lat+0.175])
            gl = ax.gridlines(draw_labels=False)
            ax.add_image(img, img_quality)
            pcm1 = ax.pcolormesh(self.lon_grid, self.lat_grid, self.delR_mbmp_final_plume_only, cmap=plt.cm.coolwarm, vmin=vmin_delR, vmax=-vmin_delR, shading='auto', transform=ccrs.PlateCarree()) 

            ####  Masks
            vmin_aux_main = np.mean(self.auxiliary_bands_main[-1]) - standard_deviations*np.std(self.auxiliary_bands_main[-1])
            vmax_aux_main = np.mean(self.auxiliary_bands_main[-1]) + standard_deviations*np.std(self.auxiliary_bands_main[-1])
            vmin_aux_ref = np.mean(self.auxiliary_bands_ref[-1]) - standard_deviations*np.std(self.auxiliary_bands_ref[-1])
            vmax_aux_ref = np.mean(self.auxiliary_bands_ref[-1]) + standard_deviations*np.std(self.auxiliary_bands_ref[-1])  
            ax = fig.add_subplot(2, 4, 4, projection=img.crs) 
            mask_title = 'Masks       '
            for i in range(len(self.bands_in_aux_mask)):
                if i != 0:
                    mask_title = mask_title + '   |   '
                mask_title = mask_title + self.bands_in_aux_mask[i]
            ax.set_title(mask_title, fontsize=10)
            ax.set_extent([self.lon-0.15, self.lon+0.15, self.lat-0.15, self.lat+0.15], crs=ccrs.PlateCarree())
#             ax.set_extent([self.lon-0.225, self.lon+0.225, self.lat-0.175, self.lat+0.175])
            gl = ax.gridlines(draw_labels=False)
            if self.filter_delR_mask:
                pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.delR_mask, cmap=plt.cm.Blues, vmin=0, vmax=1., alpha=1., shading='auto', transform=ccrs.PlateCarree())
            if self.filter_auxiliary_mask:
                pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.auxiliary_mask, cmap=plt.cm.Reds, vmin=0, vmax=1., alpha=0.1, shading='auto', transform=ccrs.PlateCarree())

             ####  Visual band main date (S1) including correlation with reference date
            ax = fig.add_subplot(2, 4, 7, projection=img.crs)
            ax.set_title('Main visual  |   $R_{vis}$=' + str(np.around(float(self.aux_corr), decimals=2)), fontsize=10)
            ax.set_extent([self.lon-0.15, self.lon+0.15, self.lat-0.15, self.lat+0.15], crs=ccrs.PlateCarree())
#             ax.set_extent([self.lon-0.225, self.lon+0.225, self.lat-0.175, self.lat+0.175])
            pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.auxiliary_bands_main[-1], vmin=vmin_aux_main, vmax=vmax_aux_main, cmap=plt.cm.gist_gray, shading='auto', transform=ccrs.PlateCarree())
            gl = ax.gridlines(draw_labels=False)
            ax.scatter(self.lon, self.lat, s=70, c='white', marker='x', transform=ccrs.PlateCarree())

            ####  Visual band reference date (S1) with wind    
            ax = fig.add_subplot(2, 4, 8, projection=img.crs)
            ax.set_extent([self.lon-0.15, self.lon+0.15, self.lat-0.15, self.lat+0.15], crs=ccrs.PlateCarree())
#             ax.set_extent([self.lon-0.225, self.lon+0.225, self.lat-0.175, self.lat+0.175])
            pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.auxiliary_bands_ref[-1], vmin=vmin_aux_ref, vmax=vmax_aux_ref, cmap=plt.cm.gist_gray, shading='auto', transform=ccrs.PlateCarree())
            ax.scatter(self.lon, self.lat, s=70, c='white', marker='x', transform=ccrs.PlateCarree())
            if type(self.ref_date) != str:
                ax.set_title('Ref. visual  |  ' + str(np.around(self.wind_ref, decimals=1)) + 'm/s', fontsize=10)
                plt.quiver(lon_wind_1, lat_wind_1, self.u10_ref_norm, self.v10_ref_norm, scale = 10, color='black', transform=ccrs.PlateCarree())
            else:
                ax.set_title('Ref. visual', fontsize=10)
            gl = ax.gridlines(draw_labels=False)

            fig.tight_layout()

            if save:          
                save_name = './StandardImages/Sentinel_3/' + self.plume_name + '_' + self.main_date.strftime('%Y-%m-%d_%H%M') + '.png'
                fig.savefig(save_name, bbox_inches='tight')
            
            if show:
                plt.show()
            else:
                plt.close()

            if save:   
                S3_save_dict = {'plume_name': self.plume_name,
                                'main_date': self.main_date,
                                'ref_date': self.ref_date,
                                'cen_lon': self.lon,
                                'cen_lat': self.lat,
                                'lon_grid': self.lon_grid,
                                'lat_grid': self.lat_grid,
                                'scene_width_height': 2*self.delta_lonlat,
                                'delR_original': self.delR_mbmp_original,
                                'delR_destriped': self.delR_mbmp_destriped,
                                'delR_filtered': self.delR_mbmp_final,
                                'omega': self.omega,
                                'main_swir2': self.swir2_main,
                                'main_swir1': self.swir1_main,
                                'ref_swir2': self.swir2_ref,
                                'ref_swir1': self.swir1_ref,
                                'mask': self.mask_to_save,
                                'sza_main': self.sza_main,
                                'vza_main': self.vza_main,
                                'sza_ref': self.sza_ref,
                                'vza_ref': self.vza_ref,
                                'wind_speed_main': self.wind_main,
                                'wind_direction_main': self.wind_dir_main,
                                'wind_speed_ref': self.wind_ref,
                                'wind_direction': self.wind_dir_ref
                               }
                
                with open('./Pickle_files/Sentinel_3/' + self.plume_name + '_' + self.main_date.strftime('%Y-%m-%d_%H%M') + '.pickle', 'wb') as handle:
                    pickle.dump(S3_save_dict, handle)
            
        else:
            print()
            print('Invalid data for Standard Image S3, showing main date data: ')
            
            fig = plt.figure(figsize=(14, 7))
            fig.suptitle(str(self.main_date.date()) + '  |  No valid reference date')
            ax = fig.add_subplot(1, 2, 1, projection=ccrs.PlateCarree())
            ax.set_title('SWIR-2/SWIR-1 ratio', loc='left',fontsize=10)
            ax.set_extent(self.grid_of_interest)
            ax.coastlines()
            gl = ax.gridlines(draw_labels=False)
            pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.swir2_main/self.swir1_main, cmap=plt.cm.coolwarm, shading='auto', transform=ccrs.PlateCarree())
            ax = fig.add_subplot(1, 2, 2, projection=ccrs.PlateCarree())
            ax.set_title('Visual band (S1)', loc='left',fontsize=10)
            ax.set_extent(self.grid_of_interest)
            ax.coastlines()
            gl = ax.gridlines(draw_labels=False)
            pcm = ax.pcolormesh(self.lon_grid, self.lat_grid, self.auxiliary_bands_main[-1], cmap=plt.cm.coolwarm, shading='auto', transform=ccrs.PlateCarree())
            plt.show()
            
    def create_standard_image_tropomi_and_save_data(self,
                                                    delta_lonlat_tropomi=2.0,
                                                    show=False,
                                                    destripe=False,
                                                    save=False
                                                    ):
        
        if self.verbose:
            print('Retrieving Tropomi data ...')
    
        makeCordbox(self.plume_name, self.lon, self.lat, delta_lonlat_tropomi)        
        makeOrbitWisePlotsCartopy(self.plume_name, self.main_date, end_day=self.main_date, aot_alb=True, destripe=destripe, show=show, save=save)        
  

    def create_standard_image_S2_and_save_data(self,
                                               area=4,
                                               show=False,
                                               search_ref_date=True,
                                               multi_refs=True,
                                               save=False
                                               ):
        
        if self.verbose:
            print('Retrieving Sentinel 2 data ...')
        
        S2_object = ms.ch4ret(self.main_date, self.main_date-td(days=20), 'MBMP', self.lat, self.lon, area=area, case=self.plume_name, TOAsat= 'S2', verbose = self.verbose )
        S2_object.delRcalc(self.main_date-td(days=20), plot_option=True, search_ref_date=search_ref_date, multi_refs=multi_refs)
        fig = S2_object.fig_delR
        S2_sza_main = S2_object.sza1
        S2_vza_main = S2_object.vza1
        
        if not multi_refs:
            S2_sza_ref = S2_object.sza2
            S2_vza_ref = S2_object.vza2
        else:
            S2_sza_ref = 'Multi'
            S2_vza_ref = 'Multi'
            
        delR_S2 = S2_object.delR
        min_from_median = np.nanmedian(delR_S2) - delR_S2.min()
        
        if self.verbose:
            print('delR minimum - median:    ', min_from_median)

        # Convert delR (radiance difference) to Omega (methane concentration)
        S2_omega= self.fullMBMP2Omega_(S2_object.delR, 'S2', S2_sza_main, vza= S2_vza_main)
      
        if save:           
            save_name = './StandardImages/Sentinel_2/' + self.plume_name + '_' + S2_object.main_date.strftime('%Y-%m-%d_%H%M') + '.png'
            fig.savefig(save_name, bbox_inches='tight')
            
        if show:
            plt.show()
        else:
            plt.close()
 
        if save:
            S2_save_dict = {'plume_name': self.plume_name,
                            'main_date': S2_object.main_date,
                            'ref_date': S2_object.ref_date,
                            'cen_lon': S2_object.lng,
                            'cen_lat': S2_object.lat,
                            'area': area,
                            'delR': S2_object.delR,
                            'lonlat_bounds': S2_object.getRegion(),
                            'omega': S2_omega,
                            'main_swir2': S2_object.R12,
                            'main_swir1': S2_object.R11,
                            'ref_swir2': S2_object.R12old,
                            'ref_swir1': S2_object.R11old,
                            'sza_main': S2_sza_main,
                            'vza_main': S2_vza_main,
                            'sza_ref': S2_sza_ref,
                            'vza_ref': S2_vza_ref,
                            'wind_speed': S2_object.wind_speed,
                            'wind_direction': S2_object.wind_dir
                           }                         

            with open('./Pickle_files/Sentinel_2/' + self.plume_name + '_' + S2_object.main_date.strftime('%Y-%m-%d_%H%M') + '.pickle', 'wb') as handle:
                pickle.dump(S2_save_dict, handle)   
                    
                    
    def create_standard_images_and_save_data_to_pickle_files(self,
                                                             S3_img_quality=10,
                                                             S3_show=False,
                                                             retrieve_tropomi=True,
                                                             tropomi_delta_lonlat=2.0,
                                                             tropomi_show=False,
                                                             tropomi_destripe=False,
                                                             retrieve_S2=True,
                                                             S2_area=4,
                                                             S2_show=False,
                                                             S2_search_ref_date=True,
                                                             S2_multi_refs=True,
                                                             save_results=False
                                                            ):
        
        self.create_standard_image_and_save_data(img_quality = S3_img_quality,
                                                 show = S3_show,
                                                 save=save_results)
        
        if retrieve_tropomi:
            self.create_standard_image_tropomi_and_save_data(delta_lonlat_tropomi=tropomi_delta_lonlat,
                                                             show=tropomi_show,
                                                             destripe=tropomi_destripe,
                                                             save=save_results
                                                             )
    
        if retrieve_S2:
            self.create_standard_image_S2_and_save_data(area=S2_area,
                                                        show=S2_show,
                                                        search_ref_date=S2_search_ref_date,
                                                        multi_refs=S2_multi_refs,
                                                        save=save_results
                                                        )
            