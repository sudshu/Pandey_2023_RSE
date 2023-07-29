# Sentinel-3 Methane Detection

Reference:
Pandey, S., van Nistelrooij, M., Maasakkers, J. D., Sutar, P., Houweling, S., Varon, D. J., Tol, P., Gains, D., Worden, J., and Aben, I.: Daily detection and quantification of methane leaks using Sentinel-3: a tiered satellite observation approach with Sentinel-2 and Sentinel-5p, Remote Sens. Environ., 296, 113716, https://doi.org/10.1016/j.rse.2023.113716, 2022. 



### Utilities and Sentinel-3 Class
Sentinel3_Utilities.py contains some standalone functions that are used throughout the processing, in several files, as discussed below.
Sentinel3Class.py contains a class object that is used to analyze the Sentinel-3 observations. This class contains several functions required for the analysis, that must be executed in order. Some of these functions have optional input parameters to finetune the process on a case-by-case basis.
Code_IME_CSF contains code for the quantification of methane emissions according to the IME and CSF methods, but this was not implemented yet. test_srf_210104_delr_to_omega.pkl contains data used for the conversion of DelR_MBMP values to methane concentrations (omega) as used in Sentinel3Class.py.

## Procedure

### Downloading
Sentinel-3 SLSTR observations are downloaded from https://scihub.copernicus.eu/dhus/#/home through the Sentinelsat API. Two scripts are available to download the data. Both scripts work very similarly, with one difference being that one of them downloads observations of dates that are provided within an CSV input file, whereas the other downloads observations in a time range. These scripts are called Download_dates_from_CSV.py and Download_single_date_range.py. Additional inputs for these scripts are given within them and are:
- head_folder, data_folder: names of directories where to store the downloads
- username, password: from https://scihub.copernicus.eu/dhus/#/home account as required
by the API
- longitude, latitude: geographical location to consider
- platformname, producttype: satellite and data package that is to be downloaded
The dates that are to be considered are either given by:
- dates_file: the filename of the CSV file containing the dates. This CSV file should have three
columns, one called Days, one called Months and one Called Years. Each row is then a date
that should be downloaded
- date_start, day_range: defining the begin and range of dates that should be downloaded.
Both scripts include a day-time filter, so they only download valid observations taken when light is present. In addition, after downloading, the scripts unzip the downloaded data files so that they can be used immediately. Both the main and reference observations can be stored in the same folder. Of the downloaded packages, not al files are required for methane detection. The ones not needed, can be deleted (if desired because of disk space condirations for example) with Remove_useless.ipynb, which only takes the directory containing the files as input.
### Interpolating
Having downloaded the data, an interpolation, or re-gridding, is required to match the pixels of main and reference observations. This is done for only a portion of the complete scene that was downloaded. This is one in a separate script and the results are saved, such that this time-consuming step is only required once. The script Interpolate.py takes the following inputs:
- cen_lon, cen_lat: central longitude and latitude of the scene that is to be considered.
- delta_lonlat: half of the width of the scene that is to be considered for interpolation
- pixel_size: size of the new pixels, default is 0.001, which was found to work well.
- case_name: name of the location, can be useful for filenaming for example
- original_file_dir: directory where the original files can be found
- interpolated_dir: directory where the interpolation results are to be stored
This interpolation can be done multiple times for one case with different centers or scene widths (delta_lonlat) if the exact source location is unknown or further zooming is required.
### Inputs
The next step is defining inputs for the processing of the data, which is done in Inputs.py. This file contains several sections, each with inputs for one specific case. First, a read me part is present, which tells you how the main dates can be defined. Then, for each case, the following inputs are combined into an input dictionary:
- plume_name: name of the location, useful for filenaming for example
- cen_lon, cen_lat: central longitude and latitude of the scene that is to be considered.
- delta_lonlat: half of the width of the scene that is to be considered for interpolation
- pixel_size: size of the new pixels, default is 0.001, which was found to work well.
- all_data_dir: directory containg the interpolated data, both main and reference
observations
- specific_main_dates: a list with datetime objects defining with dates are to be considered as
main dates, that is potentially containing methane plumes. There are three ways to define the specific_main_dates, as explained at the beginning of Inputs.py.
### Processing
Process.py and Process_S3_only.py are subsequently used to analyze the Sentinel-3 observations. These scripts make extensive use of the Sentinel3Class.py. The difference is that Process.py analyzes, besides Sentinel-3, the corresponding TROPOMI and Sentinel-2 observations as well. These scripts can save the relevant results as Pickle files, as well as the developed standard images to assess the outcome (also for TROPOMI and Sentinel-2).
These scripts require the following inputs:
- min_vis_corr: minimum correlation with main observation in the visible band (called S1) to
consider a reference observation
- min_delR: minimum value of DelR_MBMP to consider the results valid
- number_of_references_original: the number of best reference observations that should be
used to construct the final reference observations by taking the pixel-wise median. If not
enough valid once are present, less will be used.
- reference_day_range: range of days ahead and before the main date that are considered
valid reference observations.
Optional: (defaults are present in Sentinel3Class.py)
- destripe: whether to apply the destriping procedure
- filter_delR_mask: whether to apply the structural similarity mask between DelR_MBSP to
identify regions of interest (potential methane plumes)

- delR_mask_limit: maximum DelR_MBSP value to consider in that filter, as to remove already signals that cannot be due to methane
- filter_auxiliary_mask: whether to apply the structural similarity mask between auxiliary bands to identify noise regions
- correlation_minimum: minimum auxiliary band correlation to consider it in the mask
- auxiliary_bands_to_consider: bands to consider in the auxiliary mask
- filter_top_xperc: whether to include a top x percentile filter
- top_xperc: the value for x in the filter
- filter_median: whether to include a median filter
- median_filter_size: size of median filter Result variables:
- SHOW_STANDARD_TROPOMI: whether to show the TROPOMI standard image
- SHOW_STANDARD_S3: whether to show the Sentinel 3 standard image
- SHOW_STANDARD_S2: whether to show the Sentinel 2 standard image
- save_results: whether to save the standard images and Pickle files
- img_quality_S3: quality of background in S3 standard image
- destripe_tropomi: whether to destripe the Tropomi image
- delta_lonlat_tropomi: half width of tropomi scene
- search_ref_date_S2: whether to do a correlation based search for Sentinel 2 reference date
- multi_refs_S2: whether to use the median of multiple reference observations for Sentinel 2
- area_S2: scene width of Sentinel 2
- VERBOSE: whether to print information during the processing
- input_dict: dictionary containing the cases with their corresponding case specific inputs as
defined in Inputs.py
