
import io, glob, os
import pickle

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

from urllib.request import urlopen, Request
from PIL import Image
#import cv2 as cv

from pylab import *
import folium
from folium.plugins import MousePosition, MiniMap
from folium import IFrame, plugins
import base64

from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import rgb2hex, Normalize, LogNorm
from matplotlib.patches import Polygon as Poly

from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import GoogleTiles
from cartopy import geodesic
from shapely.geometry import LineString, Polygon, mapping

import scipy.spatial.qhull as qhull
import netCDF4 as nc


################################################################################


def read_full_slstr_data(filepath):
    
    """
    Read all required Sentinel 3 SLSTR data from a certain file, it takes the cn files, which is a combination of both stripe A and B nadir view.
    
    Inputs:
    - filepath: str
    
    Returns:
    - lon: np.ndarray
    - lat: np.ndarray
    - s6: np.ndarray
    - s5: np.ndarray
    - s4: np.ndarray
    - s3: np.ndarray
    - s2: np.ndarray
    - s1: np.ndarray
    
    - landscape: np.ndarray
    
    """
    
    
    data_types = ['/geodetic_an.nc', '/S6_radiance_an.nc', '/S5_radiance_an.nc',
                  '/S4_radiance_an.nc', '/S3_radiance_an.nc', '/S2_radiance_an.nc',
                  '/S1_radiance_an.nc', '/flags_an.nc', '/met_tx.nc', 'geometry_tn.nc',
                  '/geodetic_tx.nc'
                 ]
    
    geodetic = nc.Dataset(filepath + '/' + data_types[0])
    s6 = nc.Dataset(filepath + '/' + data_types[1])
    s5 = nc.Dataset(filepath + '/' + data_types[2])
    s4 = nc.Dataset(filepath + '/' + data_types[3])
    s3 = nc.Dataset(filepath + '/' + data_types[4])
    s2 = nc.Dataset(filepath + '/' + data_types[5])
    s1 = nc.Dataset(filepath + '/' + data_types[6])

    flags = nc.Dataset(filepath + '/' + data_types[7])
    mettx = nc.Dataset(filepath + '/' + data_types[8])
    geometry = nc.Dataset(filepath + '/' + data_types[9])
    geodetic_tx = nc.Dataset(filepath + '/' + data_types[10])
        
    lon = geodetic['longitude_an'][:]
    lat = geodetic['latitude_an'][:]
    s6 = s6['S6_radiance_an'][:]
    s5 = s5['S5_radiance_an'][:]
    s4 = s4['S4_radiance_an'][:]
    s3 = s3['S3_radiance_an'][:]
    s2 = s2['S2_radiance_an'][:]
    s1 = s1['S1_radiance_an'][:]
    
    cloud = flags['cloud_an'][:]
    landscape = flags['confidence_an'][:]
    
    cloud_frac = mettx['cloud_fraction_tx'][:]
    
    sza = geometry['solar_zenith_tn'][:]
    vza = geometry['sat_zenith_tn'][:]
    longitude_tx = geodetic_tx['longitude_tx'][:]
    latitude_tx = geodetic_tx['latitude_tx'][:]

    return lon, lat, s6, s5, s4, s3, s2, s1, landscape, sza, vza, longitude_tx, latitude_tx


def read_full_interpolated_data(filepath):

    """

    Read the interpolated S6, S5, S4, S3, S2 and S1 band, water mask data and Solar Zenith Angle and Viewing Zenith Angle

    """

    s6 = np.load(filepath + '/S6_interpolated.npy')
    s5 = np.load(filepath + '/S5_interpolated.npy')
    s4 = np.load(filepath + '/S4_interpolated.npy')
    s3 = np.load(filepath + '/S3_interpolated.npy')
    s2 = np.load(filepath + '/S2_interpolated.npy')
    s1 = np.load(filepath + '/S1_interpolated.npy')
    water_mask = np.load(filepath + '/water_interpolated.npy')

    sza = np.load(filepath + '/SZA_average.npy')
    vza = np.load(filepath + '/VZA_average.npy')

    return s6, s5, s4, s3, s2, s1, water_mask, sza, vza


def get_data_on_grid_of_interest(lon,
                                 lat,
                                 z_list,
                                 bounds):
    
    """
    Reduce data arrays to correspond to a certain region, restricted by longitude and 
    latitude boundaries
    
    Inputs:
    - lon: np.ndarray
    - lat: np.ndarray
    - z_list: np.ndarray
    - bounds: list with 4 floats
    
    Returns:
    - lon: np.ndarray
    - lat: np.ndarray
    - z_list: np.ndarray
    
    """
    
    indices = np.array(np.where(np.logical_and(np.logical_and(lon>=bounds[0], lon<=bounds[1]), np.logical_and(lat>=bounds[2], lat<=bounds[3]) ) ) )
        
    indices1 = np.unique(indices[0])
    indices2 = np.unique(indices[1])
    grid = np.meshgrid(indices1, indices2)
    
    lon = lon[grid]
    lat = lat[grid]
    for i in range(len(z_list)):
        z_list[i] = z_list[i][grid]
    return lon, lat, z_list


def get_interpolation_weights(lon_old, lat_old, lon_new, lat_new, d=2):
    
    """
    Get weights and vertices to interpolate from one grid to another
    
    Inputs:
    - lon_old: np.ndarray
    - lat_old: np.ndarray
    - lon_new: np.ndarray
    - lat_new: np.ndarray
    - d: int, optional
    
    Returns:
    - vertices
    - weights
    
    """
    
    if len(lon_old.shape) == 2:
        lonlat_old = np.zeros([lon_old.shape[0]*lon_old.shape[1],2])
    else:
        lonlat_old = np.zeros([lon_old.shape[0], 2])
    lonlat_old[:,0] = lon_old.flatten()
    lonlat_old[:,1] = lat_old.flatten()
    lonlat_new = np.zeros([lon_new.shape[0]*lon_new.shape[1],2])
    lonlat_new[:,0] = lon_new.flatten()
    lonlat_new[:,1] = lat_new.flatten()
    xyz = lonlat_old
    uvw = lonlat_new
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


def interpolate(values, vtx, wts, fill_value=0):
    
    """
    Perform interpolation on array with new vertices and weights
    
    Inputs:
    - values: np.ndarray
    - vtx: vertices
    - wts: weights
    
    Returns:
    - ret: np.ndarray
    
    """
    
    values=values.flatten()
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < -1e-7, axis=1)] = fill_value
    return ret


def perform_interpolation(x_old, y_old, x_new, y_new, z_list):
    
    """
    Perform interpolation, with get_inteporlation_weights() and interpolate(), on array from old to new grid
    
    Inputs:
    - x_old: np.ndarray
    - y_old: np.ndarray
    - x_new: np.ndarray
    - y_new: np.ndarray
    - z_list: np.ndarray
    
    Returns:
    - z_list: np.ndarray
    
    """
    
    vtx, wts = get_interpolation_weights(x_old, y_old, x_new, y_new)
 
    for i in range(len(z_list)):
        z_list[i] = np.reshape(interpolate(z_list[i], vtx, wts), (len(x_new), y_new[0].size))
        
    return z_list


def get_inland_water_mask(landscape_original):
    
    """
    Convert landscape values into an inland water mask
    
    Inputs:
    - landscape_original: np.ndarray
    
    Returns:
    - water_mask: np.ndarray
    
    """
    
    original_shape = landscape_original.shape
    unique_ids = np.unique(landscape_original)
    water_ids = []
    for value in unique_ids:
        flags_present = []
        rest = value  
        while rest > 2:
            power = (np.floor(np.log2(rest)))
            flag = 2**power
            flags_present.append(flag)
            rest = rest - flag  
        if 16.0 in flags_present:
            water_ids.append(value)        
            
    water_mask = np.ones(original_shape)
    for water_id in water_ids:
        water_mask[ np.where(landscape_original == water_id)] = 0

    return water_mask


def normalize_negative_one(img):

    """

    Normalized the array in the range [-1,1]

    Input:
    - Array

    Output:
    - Normalized array

    """

    normalized_input = (img - np.amin(img)) / np.ptp(img)
    return 2*normalized_input - 1

def normalize(img):

    """

    Normalized the array in the range [0,1]

    Input:
    - Array

    Output:
    - Normalized array

    """

    normalized_input = (img - np.amin(img)) / np.ptp(img)

    return normalized_input


def to_integer(img):

    """

    Converts array with floating point data type in the range [0,1]

    Input:
    - Array

    Output:
    - Normalized array
    Function to convert an array with a floating point data type in the range [0,1]
    to 8-bit unsigned integers in the range [0,255]
    input: np array normalized in the range [0,1]
    output: converted np array
    """


    integer_input = (img * 255).astype("uint8")
    return integer_input

    
def set_nan_inf_neg_to_zero(quantity_list):
    """
    Filter all not-a-number and (-) infinities in all arrays and set them to zero
    """    
    for i in range(len(quantity_list)):
        quantity_list[i][np.isinf(quantity_list[i])] = 0.
        quantity_list[i][np.isnan(quantity_list[i])] = 0.
        quantity_list[i][np.where(quantity_list[i] < 0)] = 0.

    return quantity_list


def create_zoom_images(plume_name,
                       cen_lat,
                       cen_lon,
                       main_date
                      ):
    
        S3_files = glob.glob('./Pickle_files/Sentinel_3/' + plume_name + '_' + main_date.strftime('%Y-%m-%d_%H%M') + '.pickle')

        tropomi_files = glob.glob('./Pickle_files/Tropomi/' + plume_name + '_' + str(main_date.date()) + '_*.pickle')    

        S2_files = glob.glob('./Pickle_files/Sentinel_2/' + plume_name + '_' + str(main_date.date()) + '_*.pickle')

        S3_file = open(S3_files[0], 'rb')
        S3_dict = pickle.load(S3_file)
        S3_file.close()
        
        if S2_files:
            S2_file = open(S2_files[0], 'rb')
            S2_dict = pickle.load(S2_file)
            S2_file.close()
        else:
            S2_dict = False
        
        if tropomi_files:
            for i in range(len(tropomi_files)):
                tropomi_file = open(tropomi_files[i], 'rb')
                tropomi_dict = pickle.load(tropomi_file)
                tropomi_file.close()
                
                create_one_tropomi_S3_S2_zoom(plume_name,
                                              cen_lon,
                                              cen_lat,
                                              S3_dict,
                                              tropomi_dict,
                                              S2_dict
                                             )

        else:
            tropomi_dict = False
            create_one_tropomi_S3_S2_zoom(plume_name,
                                          cen_lon,
                                          cen_lat,
                                          S3_dict,
                                          tropomi_dict,
                                          S2_dict
                                         )


def create_one_tropomi_S3_S2_zoom(plume_name,
                                  cen_lon,
                                  cen_lat,
                                  S3_dict,
                                  tropomi_dict,
                                  S2_dict
                                 ):  
    
    standard_deviations = 3
            
    if cen_lon >= 0.:
        lon_unit = '$\degree$E'
    else:
        lon_unit = '$\degree$W'
    if cen_lat >= 0.:
        lat_unit = '$\degree$N'
    else:
        lat_unit = '$\degree$S'   

    omega_S3 = S3_dict['omega']
    lon_grid_S3 = S3_dict['lon_grid']
    lat_grid_S3  = S3_dict['lat_grid']
    date_S3 = S3_dict['main_date']
    wind_speed_S3 = S3_dict['wind_speed_main']
    wind_dir_S3 = S3_dict['wind_direction_main']
    delta_lonlat_S3 = S3_dict['scene_width_height']

    vmax_omega_S3 = np.nanmean(omega_S3) + standard_deviations * np.nanstd(omega_S3)
    
    # Get 0.1 deg distance in km S3
    myGeod = geodesic.Geodesic()
    shapelyObject_horizontal = LineString([(cen_lon, cen_lon+0.1), (cen_lat, cen_lat)])
    horizontal_distance_S3 = int(myGeod.geometry_length(shapelyObject_horizontal)/1000)

    if tropomi_dict:
        omega_tropomi = tropomi_dict['xch4']
        date_tropomi = tropomi_dict['main_date']
        time_tropomi = tropomi_dict['time']
        delta_lonlat_tropomi = tropomi_dict['scene_width_height']
        lon_corners_tropomi = tropomi_dict['lon_corners']
        lat_corners_tropomi = tropomi_dict['lat_corners']
        
        vmin_omega_tropomi = np.nanmean(omega_tropomi) - standard_deviations * np.nanstd(omega_tropomi)

        # Tropomi box coordinates
        tropomi_box_left_x = cen_lon - delta_lonlat_S3/2
        tropomi_box_right_x = cen_lon + delta_lonlat_S3/2
        tropomi_box_upper_y = cen_lat + delta_lonlat_S3/2
        tropomi_box_lower_y = cen_lat - delta_lonlat_S3/2

        # Get 1 deg distance in km tropomi
        myGeod = geodesic.Geodesic()
        shapelyObject_horizontal = LineString([(cen_lon, cen_lon+1.0), (cen_lat, cen_lat)])
        horizontal_distance_trop = int(myGeod.geometry_length(shapelyObject_horizontal)/1000)
        
    if S2_dict:
        omega_S2 = S2_dict['omega']
        date_S2 = S2_dict['main_date']
        wind_speed_S2 = S2_dict['wind_speed']
        wind_dir_S2 = S2_dict['wind_direction']
        area_S2 = S2_dict['area']
        lonlat_bounds_S2 = S2_dict['lonlat_bounds']
    
        omega_S2[np.where(omega_S2 > 50)] = 0.
        vmax_omega_S2 = np.nanmean(omega_S2) + 3 * np.nanstd(omega_S2)
                
        # S3 box coordinates
        s3_box_left_x = lonlat_bounds_S2[0][0]
        s3_box_right_x = lonlat_bounds_S2[1][0]
        s3_box_upper_y = lonlat_bounds_S2[0][1]
        s3_box_lower_y = lonlat_bounds_S2[2][1]
                        
    # Create actual figure
    if tropomi_dict and S2_dict:
        fig = plt.figure(figsize=(15,4))
    else: 
        fig = plt.figure(figsize=(12,5))
        
    fig.suptitle(plume_name + '  |  ' + str(date_S3.date()) + '  |  '+ str(np.around(cen_lat, decimals=4)) + lat_unit + ' ' + str(np.around(cen_lon, decimals=4)) + lon_unit + '\n', fontsize=10)
    
    #  Tropomi plot
    if tropomi_dict:
        if S2_dict:
            ax1 = fig.add_subplot(1,4,1, projection=ccrs.PlateCarree())
        else:
            ax1 = fig.add_subplot(1,3,1, projection=ccrs.PlateCarree())
            
        ax1.set_title('TROPOMI  |  ' + str(time_tropomi)[0:5], fontsize=10)

        ax1.set_extent([cen_lon-delta_lonlat_tropomi/2, cen_lon+delta_lonlat_tropomi/2, cen_lat-delta_lonlat_tropomi/2, cen_lat+delta_lonlat_tropomi/2])
        cmap1 = plt.cm.inferno_r
        cmap1.set_bad(color='grey')
        norm1 = Normalize(vmin=np.mean(omega_tropomi)-2*np.std(omega_tropomi), vmax= np.mean(omega_tropomi)+2*np.std(omega_tropomi))
        mapper1 = ScalarMappable(norm=norm1, cmap=cmap1)
        for j in range(len(lon_corners_tropomi)):  
            color=rgb2hex(mapper1.to_rgba(omega_tropomi[j]))
            lat_long = [(lon_corners_tropomi[j][k], lat_corners_tropomi[j][k]) for k in range(4)]
            ax1.add_patch(Poly(lat_long, closed=True,  facecolor=color, edgecolor=color))
        
#         cb1 = plt.colorbar(mapper1, shrink=0.6, orientation='horizontal')
#         cb1.set_label('[ppb]')
        ax1.set_axis_off()
        distance1 = ConnectionPatch(xyA=(cen_lon-0.45*delta_lonlat_tropomi, cen_lat-0.54*delta_lonlat_tropomi), coordsA=ax1.transData, 
                               xyB=(cen_lon-0.45*delta_lonlat_tropomi+1.0, cen_lat-0.54*delta_lonlat_tropomi), coordsB=ax1.transData, 
                                  color = 'black')
        fig.add_artist(distance1)    
        ax1.text(cen_lon-0.45*delta_lonlat_tropomi+0.15, cen_lat-0.62*delta_lonlat_tropomi, str(horizontal_distance_trop) + ' km', fontsize=12, fontweight='semibold', color='black')
    
    # Sentinel 3 plot
    if tropomi_dict and S2_dict:
        ax2 = fig.add_subplot(1,4,2, projection=ccrs.PlateCarree())   
    elif tropomi_dict and not S2_dict:
        ax2 = fig.add_subplot(1,3,2, projection=ccrs.PlateCarree())   
    elif not tropomi_dict and S2_dict:
        ax2 = fig.add_subplot(1,3,1, projection=ccrs.PlateCarree())   

    ax2.set_title(str(np.around(abs(wind_speed_S3), decimals=1)) + ' m/s    |    Sentinel 3    |    ' + str(date_S3.time())[0:5], fontsize=10)

    omega_S3_plot = omega_S3.copy()
    omega_S3_plot[np.where(omega_S3<0)] = 0
    pcm2 = ax2.pcolormesh(lon_grid_S3, lat_grid_S3, omega_S3_plot, cmap=plt.cm.inferno_r, vmin=0, vmax=vmax_omega_S3, shading='auto', transform=ccrs.PlateCarree())     
#     cb2 = plt.colorbar(pcm2, shrink=0.6, orientation='horizontal')
#     cb2.set_label('[ppb]')
    ax2.set_extent([cen_lon-delta_lonlat_S3/2, cen_lon+delta_lonlat_S3/2, cen_lat-delta_lonlat_S3/2, cen_lat+delta_lonlat_S3/2])

    s3_x_min = ax2.get_xlim()[0]; s3_x_max = ax2.get_xlim()[1]; s3_y_min = ax2.get_ylim()[0]; s3_y_max = ax2.get_ylim()[1];
    s3_arrow_x = cen_lon - 0.35*delta_lonlat_S3
    s3_arrow_y = cen_lat + 0.35*delta_lonlat_S3

    ax2.quiver(s3_arrow_x, s3_arrow_y, 0.8, 0.8, angles=np.array([wind_dir_S3]), scale = 10, color='black')
    ax2.scatter(cen_lon, cen_lat, s=40, c='white', marker='x')
    ax2.set_axis_off()
    
    distance2 = ConnectionPatch(xyA=(cen_lon-0.45*delta_lonlat_S3, cen_lat-0.54*delta_lonlat_S3), coordsA=ax2.transData, 
                               xyB=(cen_lon-0.45*delta_lonlat_S3+0.1, cen_lat-0.54*delta_lonlat_S3), coordsB=ax2.transData, 
                                  color = 'black')
    fig.add_artist(distance2) 
    ax2.text(cen_lon-0.45*delta_lonlat_S3+0.025, cen_lat-0.62*delta_lonlat_S3, str(horizontal_distance_S3) + ' km', fontsize=14, fontweight='semibold', color='black')

    # Sentinel 2 plot
    if S2_dict:
        if tropomi_dict:
            ax3 = fig.add_subplot(1,4,3)
        else:
            ax3 = fig.add_subplot(1,3,2)
            
        ax3.set_title(str(np.around(abs(wind_speed_S2), decimals=1)) + ' m/s    |    Sentinel 2    |    ' + str(date_S2.time())[0:5], fontsize=10)

        omega_S2_plot = omega_S2.copy()
        omega_S2_plot[np.where(omega_S2<0)] = 0

        pcm3 = ax3.imshow(omega_S2_plot, cmap=plt.cm.inferno_r, vmin=0, vmax=vmax_omega_S2)    
#         cb3 = plt.colorbar(pcm3, shrink=0.6, orientation='horizontal')
#         cb3.set_label('[ppb]')
        
        # S2 boundary values
        s2_x_min = ax3.get_xlim()[0]; s2_x_max = ax3.get_xlim()[1]; s2_y_min = ax3.get_ylim()[0]; s2_y_max = ax3.get_ylim()[1]
        s2_x_middle = (s2_x_min + s2_x_max) / 2
        s2_y_middle = (s2_y_min + s2_y_max) / 2
        s2_arrow_x = s2_x_middle - 0.35*(s2_x_max-s2_x_min)
        s2_arrow_y = s2_y_middle + 0.35*(s2_y_max-s2_y_min)
        ax3.quiver(s2_arrow_x, s2_arrow_y, 0.8, 0.8, angles=np.array([wind_dir_S2]), scale = 10, color='black')
        ax3.scatter(s2_y_middle, s2_x_middle, s=40, c='white', marker='x')
        ax3.set_xticks([])
        ax3.set_yticks([])
        distance3 = ConnectionPatch(xyA=(s2_x_middle-0.45*(s2_x_max-s2_x_min), s2_y_middle-0.54*(s2_y_max-s2_y_min)), coordsA=ax3.transData, 
                               xyB=(s2_x_middle-0.2*(s2_x_max-s2_x_min), s2_y_middle-0.54*(s2_y_max-s2_y_min)), coordsB=ax3.transData, 
                                  color = 'black')
        fig.add_artist(distance3)
        ax3.text(s2_x_middle-0.4*(s2_x_max-s2_x_min), s2_y_middle-0.62*(s2_y_max-s2_y_min), '1 km', fontsize=14, fontweight='semibold', color='black')
    
        # Draw connecting lines and rectangles S3 and S2
        con1 = ConnectionPatch(xyA=(s3_box_right_x, s3_box_lower_y), coordsA=ax2.transData, 
                               xyB=(s2_x_min, s2_y_max), coordsB=ax3.transData, color = 'black')
        con2 = ConnectionPatch(xyA=(s3_box_right_x, s3_box_upper_y), coordsA=ax2.transData, 
                               xyB=(s2_x_min, s2_y_min), coordsB=ax3.transData, color = 'black')
        fig.add_artist(con1); fig.add_artist(con2); 
        rect2 = Rectangle((s3_box_left_x, s3_box_lower_y), s3_box_right_x-s3_box_left_x, s3_box_upper_y-s3_box_lower_y, linewidth=1, edgecolor='black', facecolor='none')
        ax2.add_patch(rect2)     
    
    
    if tropomi_dict:
        # Draw connecting lines and rectangles Tropomi and S3
        con3 = ConnectionPatch(xyA=(tropomi_box_right_x, tropomi_box_lower_y), coordsA=ax1.transData, 
                               xyB=(s3_x_min, s3_y_min), coordsB=ax2.transData, color = 'black')                    
        con4 = ConnectionPatch(xyA=(tropomi_box_right_x, tropomi_box_upper_y), coordsA=ax1.transData, 
                               xyB=(s3_x_min, s3_y_max), coordsB=ax2.transData, color = 'black')
        fig.add_artist(con3); fig.add_artist(con4)
        rect1 = Rectangle((tropomi_box_left_x, tropomi_box_lower_y), delta_lonlat_S3, delta_lonlat_S3, linewidth=1, edgecolor='black', facecolor='none')
        ax1.add_patch(rect1)
        
  
    # Google Earth image
    cimgt.QuadtreeTiles.get_image = image_spoof # reformat web request for street map spoofing
    img = cimgt.QuadtreeTiles() # spoofed, downloaded street map
    
    if tropomi_dict and S2_dict:
        ax4 = fig.add_subplot(1,4,4, projection=ccrs.PlateCarree())
        ax4.set_title('Google Earth image')
        lon_min = lonlat_bounds_S2[0][0] + (lonlat_bounds_S2[1][0]-lonlat_bounds_S2[0][0])*0.35
        lon_max = lonlat_bounds_S2[0][0] + (lonlat_bounds_S2[1][0]-lonlat_bounds_S2[0][0])*0.65
        lat_min = lonlat_bounds_S2[0][1] + (lonlat_bounds_S2[2][1]-lonlat_bounds_S2[0][1])*0.35
        lat_max = lonlat_bounds_S2[0][1] + (lonlat_bounds_S2[2][1]-lonlat_bounds_S2[0][1])*0.65
        ax4.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax4.add_image(img, 14)
        ax4.scatter(cen_lon, cen_lat, s=40, c='white', marker='x')
        
    elif S2_dict:
        ax4 = fig.add_subplot(1,3,3, projection=ccrs.PlateCarree())
        ax4.set_title('Google Earth image')
        lon_min = lonlat_bounds_S2[0][0] + (lonlat_bounds_S2[1][0]-lonlat_bounds_S2[0][0])*0.35
        lon_max = lonlat_bounds_S2[0][0] + (lonlat_bounds_S2[1][0]-lonlat_bounds_S2[0][0])*0.65
        lat_min = lonlat_bounds_S2[0][1] + (lonlat_bounds_S2[2][1]-lonlat_bounds_S2[0][1])*0.35
        lat_max = lonlat_bounds_S2[0][1] + (lonlat_bounds_S2[2][1]-lonlat_bounds_S2[0][1])*0.65
        ax4.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax4.add_image(img, 14)
        ax4.scatter(cen_lon, cen_lat, s=40, c='white', marker='x')

    elif tropomi_dict:
        ax4 = fig.add_subplot(1,3,3, projection=ccrs.PlateCarree())
        ax4.set_title('Google Earth image')
        ax4.set_extent([cen_lon-delta_lonlat_S3/4, cen_lon+delta_lonlat_S3/4, cen_lat-delta_lonlat_S3/4, cen_lat+delta_lonlat_S3/4])
        ax4.add_image(img, 14)
        ax4.scatter(cen_lon, cen_lat, s=40, c='white', marker='x')

    cases = glob.glob('./Trop_S3_S2_ZoomIns/' + plume_name + '_' + date_S3.strftime('%Y-%m-%d') + '*.png')
    
    if len(cases)==0:
        idx = ''
    else:
        idx = '_' + str(len(cases)+1)

    save_name = './Trop_S3_S2_ZoomIns/' + plume_name + '_' + date_S3.strftime('%Y-%m-%d') + idx + '.png'

    fig.savefig(save_name)

    if show:
        plt.show()
    else:
        plt.close()



def get_specific_date_files(dates, directory):
    
    """
    Obtain a list of file paths containing all files that contain specific dates within a certain directory
    
    Inputs:
    - dates: datetime.datetime objects
    - directory: str
    
    Returns:
    - new_file_list: list of str
    
    """
    new_file_list = []
    for date in dates:        
        if date.day < 10:
            if date.month < 10:
                new_files = glob.glob(directory + '/*__' + str(date.year) + '0' + str(date.month) + '0' + str(date.day) + 'T*')
            else:
                new_files = glob.glob(directory + '/*__' + str(date.year) + str(date.month) + '0' + str(date.day) + 'T*')
        else:
            if date.month < 10:
                new_files = glob.glob(directory + '/*__' + str(date.year) + '0' + str(date.month) + str(date.day) + 'T*')
            else:
                new_files = glob.glob(directory + '/*__' + str(date.year) + str(date.month) + str(date.day) + 'T*')
                
        for new_file in new_files:
            new_file_list.append(new_file)
               
    return new_file_list
    
def image_spoof(self, tile):
    
    """
    Loads background image for Standard Image Sentinel 3 to be used with Cartopy
    
    Inputs:
    - tile
    
    Returns:
    - img
    - tile-extent
    
    """
    
    '''this function reformats web requests from OSM for cartopy
    Heavily based on code by Joshua Hrisko at:
        https://makersportal.com/blog/2020/4/24/geographic-visualizations-in-python-with-cartopy'''

    url = self._image_url(tile)                # get the url of the street map API
    req = Request(url)                         # start request
    req.add_header('User-agent','Anaconda 3')  # add user agent to request
    fh = urlopen(req) 
    im_data = io.BytesIO(fh.read())            # get image
    fh.close()                                 # close url
    img = Image.open(im_data)                  # open image with PIL
    img = img.convert(self.desired_tile_form)  # set image format
    return img, self.tileextent(tile), 'lower' # reformat for cartopy


def apply_wind_rotation(wind_dir,
                        delR,
                        translate_first=False,
                        show_plots=False
                       ):
    """
    Rotates an array (delR -> delR_rotated) by specified angle (wind_dir), with option to perform translation first and plot the array
    before and after the rotation.
    
    Inputs:
    - wind_dir: float
    - delR: np.ndarray
    
    Returns:
    - delR_rotated: np.ndarray
    
    """
    
    rotation_angle_towards_north = wind_dir - 90.0
    delR_to_rotate = delR.copy()
    
    if translate_first:
        translation_range_x = translate_first[0]
        translation_range_y = translate_first[1]
                
        delR_to_rotate = cv.warpAffine(delR_to_rotate, np.float32([[1,0,translation_range_x], [0,1,translation_range_y]]), delR_to_rotate.shape)

    delR_rotated = ndimage.rotate(delR_to_rotate, rotation_angle_towards_north, reshape=False)
    
    # Remove positives
    delR_rotated[np.where(delR_rotated > 0)] = 0   
    
    if show_plots:
        vmin = np.mean(delR)-3*np.std(delR)
        vmax = np.mean(delR)+3*np.std(delR)
        
        fig = plt.figure()        
        if translate_first:
            fig.suptitle('Translated first')
        else:
            fig.suptitle('Center rotation')
            
        ax = fig.add_subplot(1,2,1)
        ax.imshow(delR_to_rotate, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax, origin='lower')
#         ax.set_xlim([145, 155])
#         ax.set_ylim([145, 155])
        ax = fig.add_subplot(1,2,2)
        ax.imshow(delR_rotated, cmap=plt.cm.coolwarm, vmin=vmin, vmax=vmax, origin='lower')
#         ax.set_xlim([145, 155])
#         ax.set_ylim([145, 155])
        
        fig.show()
    
    return delR_rotated


def calc_gamma_wind_rotation(rotated_averaged
                            ):
    """
    Calculates difference (gamma) between upwind and downwind values
    
    Inputs:
    rotated_averaged: np.ndarray
    
    Returns:
    - down_wind_averaged: float
    - upwind_averaged: float
    - gamma: float
    
    """
    
    number_of_y_pixels = len(rotated_averaged)
    number_of_x_pixels = len(rotated_averaged[0])
    
    y_upwind_end = int(np.floor(number_of_y_pixels/2))
    y_downwind_start = int(np.ceil(number_of_y_pixels/2))

    # The values below can be tuned/modified
    y_upwind_start = int(np.floor(1 * number_of_y_pixels / 5))
    y_downwind_end = int(np.ceil(4 * number_of_y_pixels / 5))

    x_start = int(np.floor(2 * number_of_x_pixels / 5))
    x_end = int(np.ceil(3 * number_of_x_pixels / 5))
    
    # Calculate gamma
    downwind_pixels = rotated_averaged[y_downwind_start:y_downwind_end, x_start:x_end]
    upwind_pixels = rotated_averaged[y_upwind_start:y_upwind_end, x_start:x_end]
    
    downwind_averaged = np.mean(downwind_pixels)
    upwind_averaged = np.mean(upwind_pixels)    
    
#     plt.figure()
#     plt.imshow(downwind_pixels, origin='lower')
#     plt.show()
#     plt.figure()
#     plt.imshow(upwind_pixels, origin='lower')
#     plt.show()
#     print(downwind_pixels.max(), downwind_pixels.min(), np.median(downwind_pixels), downwind_averaged)
#     print(upwind_pixels.max(), upwind_pixels.min(), np.median(upwind_pixels), upwind_averaged)
    
    gamma = downwind_averaged - upwind_averaged
    
    return downwind_averaged, upwind_averaged, gamma


def rotate_and_average_center(plume_name,
                                 all_wind_dirs,
                                 all_delRs,
                                 all_main_dates,
                                 show_procedure=False,
                                 verbose=False,
                                 save=False
                                ):
    """
    Rotates multiple array in their center, each with each own angle, and average the resulting arrays per pixel. Includes 
    options to plot array before and after rotation and print information
    
    Inputs:
    - plume_name: str
    - all_wind_dirs: np.ndarray of floats
    - all_delRs: np.ndarray with np.ndarrays
    - show_procedure: bool
    - vebose: bool
    - save: bool
    
    Possible output:
    - Plots, shown or saved 
    - Print information
    
    """
        
    all_rotated_center = []
    for i in range(len(all_delRs)):

        wind_dir = all_wind_dirs[i]
        original = all_delRs[i]

        rotated = apply_wind_rotation(wind_dir, original, show_plots=show_procedure)
        all_rotated_center.append(rotated)

    all_rotated_center = np.array(all_rotated_center)
    center_rotated_and_averaged = np.mean(all_rotated_center, axis=0)   

    vmin = np.mean(center_rotated_and_averaged) - 3*np.std(center_rotated_and_averaged)            

    fig = plt.figure(figsize=(7,7))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_title('Average of plumes rotated in center  |  ' + str(len(all_delRs)) + ' scenes')
    pcm = ax1.imshow(center_rotated_and_averaged, cmap=plt.cm.coolwarm, vmin=vmin, vmax=-vmin, origin='lower')
    ax1.scatter(len(center_rotated_and_averaged[0])/2,len(center_rotated_and_averaged)/2, s=100, c='black', marker='x')

    if save:
        if not os.path.isdir('./Rotated_Center/' + plume_name):
            os.makedirs('./Rotated_Center/' + plume_name)
        save_name = 'Rotated_Center/' + plume_name + '/' + str(len(all_delRS)) + '_dates.png'
        fig.savefig(save_name, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

    downwind_average, upwind_average, gamma = calc_gamma_wind_rotation(center_rotated_and_averaged)

    if verbose:
        print()
        print('###')
        print()
        print('Number of rotated images averaged: ', len(all_delRs))
        print()
        print('Downwind average: ', downwind_average)
        print('Upwind average: ', upwind_average)
        print('Center rotated gamma: ', gamma)
        
        
def rotate_and_average_NON_center(plume_name,
                                  all_wind_dirs,
                                  all_delRs,
                                  all_main_dates,
                                  x_translation=100,
                                  y_translation=0,
                                  show_procedure=False,
                                  verbose=False,
                                  save=False
                                ):
    """
    Rotates multiple array after a translation first, each with each own angle, and average the resulting arrays per pixel. Includes 
    options to plot array before and after rotation and print information
    
    Inputs:
    - plume_name: str
    - all_wind_dirs: np.ndarray of floats
    - all_delRs: np.ndarray with np.ndarrays
    - x_translation: int
    - y_translation: int
    - show_procedure: bool
    - vebose: bool
    - save: bool
    
    Possible output:
    - Plots, shown or saved 
    - Print information
    
    """
    
    all_rotated_non_center = []
    
    for i in range(len(all_delRs)):
        wind_dir = all_wind_dirs[i]
        original = all_delRs[i]

        rotated = apply_wind_rotation(wind_dir, original, translate_first=[x_translation, y_translation], show_plots=show_procedure)
        all_rotated_non_center.append(rotated)
        
    all_rotated_non_center = np.array(all_rotated_non_center)
    non_center_rotated_and_averaged = np.mean(all_rotated_non_center, axis=0)   

    vmin = np.mean(non_center_rotated_and_averaged) - 3*np.std(non_center_rotated_and_averaged) 
    
    fig = plt.figure(figsize=(7,7))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_title('Average rotated plumes after (' + str(x_translation) + ',' + str(y_translation) + ') translation  |  ' + str(len(all_delRs)) + ' scenes')
    pcm = ax1.imshow(non_center_rotated_and_averaged, cmap=plt.cm.coolwarm, vmin=vmin, vmax=-vmin, origin='lower')
    ax1.scatter(len(non_center_rotated_and_averaged[0])/2,len(non_center_rotated_and_averaged)/2, s=100, c='black', marker='x')

    if save:
        if not os.path.isdir('./Rotated_NON_Center/' + plume_name):
            os.makedirs('./Rotated_NON_Center/' + plume_name)
        save_name = 'Rotated_NON_Center/' + plume_name + '/' + str(len(all_delRS)) + '_dates.png'
        fig.savefig(save_name, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        
    downwind_average, upwind_average, gamma = calc_gamma_wind_rotation(non_center_rotated_and_averaged)

    if verbose:
        print()
        print('###')
        print()
        print('Downwind average: ', downwind_average)
        print('Upwind average: ', upwind_average)
        print('NON Center rotated gamma: ', gamma)
    
    
def rotate_at_multi_pixels(plume_name,
                           all_wind_dirs,
                           all_delRs,
                           all_main_dates,
                           show_procedure=False,
                           verbose=False,
                           save=False
                          ):
    
    """
    Rotates multiple arrays at mulitple locations (regularly space in the central part of the array, each with each own angle, and averages the resulting arrays per pixel per rotation location. Includes 
    options to plot array before and after rotation and print information
    
    Inputs:
    - plume_name: str
    - all_wind_dirs: np.ndarray of floats
    - all_delRs: np.ndarray with np.ndarrays
    - all_main_dates: np.ndarray with datetime.datetime objects
    - show_procedure: bool
    - vebose: bool
    - save: bool
    
    Possible output:
    - Plots, shown or saved 
    - Print information
    
    """
    
    number_of_y_pixels = len(all_delRs[0])
    number_of_x_pixels = len(all_delRs[0,0])
    
    x_middle = number_of_x_pixels / 2
    y_middle = number_of_y_pixels / 2
    
    downwind_average_grid = np.zeros((number_of_y_pixels, number_of_x_pixels))
    upwind_average_grid = np.zeros((number_of_y_pixels, number_of_x_pixels))
    gamma_grid = np.zeros((number_of_y_pixels, number_of_x_pixels))
    
    # Loop over all pixels to rotate and average all images - start lower left, going right first and then up
    for i in np.arange(100, 201, 10): # Consider  first row (y)
        i=int(i)
        for j in np.arange(100, 201, 10): # Consider first 5 pixels (x)
            j=int(j)
            translation_range_x = x_middle - j
            translation_range_y = y_middle - i
            
            all_rotated = []
            
            # Loop over all originals / main dates
            for k in range(len(all_delRs)):
                original = all_delRs[k]
                wind_dir = all_wind_dirs[k]
                
                rotated = apply_wind_rotation(wind_dir, original, translate_first=[translation_range_x, translation_range_y], show_plots=show_procedure)
                                                
                all_rotated.append(rotated)
                            
            all_rotated = np.array(all_rotated)
            rotated_and_averaged = np.mean(all_rotated, axis=0)   
            
            downwind_average, upwind_average, gamma = calc_gamma_wind_rotation(rotated_and_averaged)
            
            downwind_average_grid[i-3:i+4, j-3:j+4] = downwind_average
            upwind_average_grid[i-3:i+4, j-3:j+4] = upwind_average
            gamma_grid[i-3:i+4, j-3:j+4] = gamma
    
    if verbose:
        print()
        print('###')
        print()
        print('Gamma grid min and max: ', gamma_grid.min(), gamma_grid.max())

    gamma_grid[np.where(gamma_grid == 0)] = 'nan'
    
    vmin = np.nanmean(gamma_grid) - 3*np.nanstd(gamma_grid)

    fig = plt.figure()
    fig.suptitle('Gamma grid  |  ' + str(len(all_delRs)) + ' scenes')
    ax = fig.add_subplot(1,1,1)
    img = ax.imshow(gamma_grid, cmap=plt.cm.coolwarm, vmin=vmin, vmax=-vmin, origin='lower')
    cbar = fig.colorbar(img, ax=ax)
    cbar.minorticks_on()
    
    if save:
        if not os.path.isdir('./Gamma_Grids/' + plume_name):
            os.makedirs('./Gamma_Grids/' + plume_name)
        save_name = 'Gamma_Grids/' + plume_name + '/' + str(len(all_delRS)) + '_dates.png'
        fig.savefig(save_name, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

        
def normalize(img):
        
        """
        
        Normalized the array in the range [0,1]
        
        Input:
        - Array
        
        Output:
        - Normalized array
        
        """
        
        normalized_input = (img - np.amin(img)) / np.ptp(img)

        return normalized_input
        
def create_folium_map(all_plume_names,
                      all_cen_lons,
                      all_cen_lats,
                      all_main_dates,
                      fig_folium_name
                     ):
    
    """
    Creates and saves a Folium Map with Tropomi, Sentinel 3 and Sentinel 2 data for multiple locations based on the saved pickle files
    
    Inputs:
    - plume_names: list of str
    - cen_lons: list of floats
    - cen_lats: list of floats
    - main_dates_S3: list of lists with datetime.datetime objects

    Returns:
    - folium map
    
    """
    
    
    # Initialize map, minimap, LatLon popup, marker layout
    tile_esri= 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
#     m= folium.Map(location=[cen_lats[0], cen_lons[0]], zoom_start=6, height=800, tiles=tile_esri,attr= 'Esri')
    m=folium.Map(height=800, tiles=tile_esri, attr='Esri')
    html = '<img src="data:image/png;base64,{}">'.format
    icon = folium.Icon(color="red", icon="ok")
    m.add_child(folium.LatLngPopup())
    MousePosition().add_to(m)
    minimap = MiniMap() ; m.add_child(minimap)     
    
    
    # Consider each location
    for i in range(len(all_plume_names)):
        plume_name = all_plume_names[i]
        print(plume_name + ' ... ')
        cen_lon = all_cen_lons[i]
        cen_lat = all_cen_lats[i]
        all_case_main_dates = all_main_dates[i]

        for j in range(len(all_case_main_dates)):
            main_date = all_case_main_dates[j]
            print(main_date)
            S3_files = glob.glob('./Pickle_files/Sentinel_3/' + plume_name + '_' + main_date.strftime('%Y-%m-%d_%H%M') + '.pickle')
            tropomi_files = glob.glob('./Pickle_files/Tropomi/' + plume_name + '_' + str(main_date.date()) + '_*.pickle')    
            S2_files = glob.glob('./Pickle_files/Sentinel_2/' + plume_name + '_' + str(main_date.date()) + '_*.pickle')
            
            ## TROPOMI
            tropomi_dicts = []
            if tropomi_files:
                if j==0:
                    for t in range(len(tropomi_files)):
                        tropomi_file = open(tropomi_files[t], 'rb')
                        tropomi_dict = pickle.load(tropomi_file)
                        tropomi_file.close()
                        tropomi_dicts.append(tropomi_dict)

                elif main_date.date() != all_case_main_dates[j-1].date():
                    for t in range(len(tropomi_files)):
                        tropomi_file = open(tropomi_files[t], 'rb')
                        tropomi_dict = pickle.load(tropomi_file)
                        tropomi_file.close()
                        tropomi_dicts.append(tropomi_dict)
                else:
                    tropomi_dicts = False
            else: 
                tropomi_dicts = False
                
            if tropomi_dicts:
                for t in range(len(tropomi_dicts)):
                    tropomi_dict = tropomi_dicts[t]

                    xch4 = normalize(tropomi_dict['xch4'])
                    lon_corners = tropomi_dict['lon_corners']
                    lat_corners = tropomi_dict['lat_corners']
                    delta_lonlat_tropomi = tropomi_dict['scene_width_height'] / 2

                    trop_fig = plt.figure()
                    ax = trop_fig.add_subplot(111, projection=ccrs.PlateCarree())
                    norm3 = Normalize(vmin=np.mean(xch4)-2*np.std(xch4), vmax= np.mean(xch4)+2*np.std(xch4))
                    mapper3 = ScalarMappable(norm=norm3, cmap=cm.jet)
                    ax.set_extent([cen_lon-delta_lonlat_tropomi, cen_lon+delta_lonlat_tropomi, cen_lat-delta_lonlat_tropomi, cen_lat+delta_lonlat_tropomi])
                    for s in range(len(lon_corners)):  
                        color=rgb2hex(mapper3.to_rgba(xch4[s]))
                        lat_long = [(lon_corners[s][r], lat_corners[s][r]) for r in range(4)]
                        CF = ax.add_patch(Poly(lat_long, closed=True,  facecolor=color, edgecolor=color))
                    trop_fig.savefig('./FoliumMaps/Tropomi_Temp' + str(t) + '.png', bbox_inches='tight')
                    plt.close()
                        
            ## Sentintel-3
            S3_file = open(S3_files[0], 'rb')
            S3_dict = pickle.load(S3_file)
            S3_file.close()

            try:
                delR_S3 = S3_dict['delR']
            except:
                delR_S3 = S3_dict['delR_original']

            delR_S3[np.where(delR_S3>0)] = 0
            vmin_S3 = np.mean(delR_S3) - 3*np.std(delR_S3)
            delta_lonlat_S3 = S3_dict['scene_width_height'] / 2

            s3_fig = plt.figure()
            ax=s3_fig.add_subplot(111, projection=ccrs.PlateCarree())
            ax.set_extent([cen_lon-delta_lonlat_S3, cen_lon+delta_lonlat_S3, cen_lat-delta_lonlat_S3, cen_lat+delta_lonlat_S3])
            ax.pcolormesh(S3_dict['lon_grid'], S3_dict['lat_grid'], delR_S3, cmap=cm.coolwarm, vmin=vmin_S3, vmax=-vmin_S3)
            s3_fig.savefig('./FoliumMaps/S3_Temp.png', bbox_inches='tight')
            plt.close()
            
            ## Sentinel-2
            if S2_files:
                if j==0: 
                    S2_file = open(S2_files[0], 'rb')
                    S2_dict = pickle.load(S2_file)
                    S2_file.close()
                elif main_date.date() != all_case_main_dates[j-1].date():
                    S2_file = open(S2_files[0], 'rb')
                    S2_dict = pickle.load(S2_file)
                    S2_file.close()
                else: 
                    S2_dict=False
            else:
                S2_dict = False
            
            if S2_dict:
                delR_S2 = S2_dict['delR']
                lonlat_bounds_S2 = S2_dict['lonlat_bounds']
                vmin_S2 = np.mean(delR_S2) - 3*np.std(delR_S2)
                
                s2_fig = plt.figure()
                ax=s2_fig.add_subplot(111)
                ax.imshow(delR_S2, cmap=cm.inferno, origin='upper')
                ax.set_xticks([])
                ax.set_yticks([])
                s2_fig.savefig('./FoliumMaps/S2_Temp.png', bbox_inches='tight')
                plt.close()
    
            ## Load into folium map
            ## Toggle on first date, toggle off the rest (per location)
            if j == 0:
                if tropomi_dicts:
                    for t in range(len(tropomi_dicts)):
                        tropomi_dict = tropomi_dicts[t]
                        tropomi_time = str(tropomi_dict['time'])[:5]
                        img_tropomi = folium.raster_layers.ImageOverlay(name=plume_name + '  |  ' + str(main_date.date()) + ' ' + tropomi_time + '  |  TROPOMI (S5p)',
                                                                        image='./FoliumMaps/Tropomi_Temp' + str(t) + '.png',
                                                                        bounds=[[cen_lat-delta_lonlat_tropomi, cen_lon-delta_lonlat_tropomi], [cen_lat+delta_lonlat_tropomi, cen_lon+delta_lonlat_tropomi]],
                                                                        origin='lower'
                                                                       ).add_to(m)
                
                img = folium.raster_layers.ImageOverlay(name=plume_name + '  |  ' + str(main_date.date()) + ' ' + str(main_date.time())[0:5] + '  |  SENTINEL-3',
                                                        image='./FoliumMaps/S3_Temp.png',
                                                        bounds=[[cen_lat-delta_lonlat_S3, cen_lon-delta_lonlat_S3], [cen_lat+delta_lonlat_S3, cen_lon+delta_lonlat_S3]],
                                                        origin='lower'
                                                        ).add_to(m)


                if S2_dict: 
                    img = folium.raster_layers.ImageOverlay(name=plume_name + '  |  ' + str(main_date.date()) + ' ' + str(S2_dict['main_date'].time())[0:5] + '  |  SENTINEL-2',
                                                            image='./FoliumMaps/S2_Temp.png',
                                                            bounds=[[lonlat_bounds_S2[0][1], lonlat_bounds_S2[0][0]], [lonlat_bounds_S2[2][1], lonlat_bounds_S2[1][0]]],
                                                            colormap=cm.inferno,
                                                            origin='lower'
                                                           ).add_to(m)
            else:
                if tropomi_dicts:
                    for t in range(len(tropomi_dicts)):
                        tropomi_dict = tropomi_dicts[t]
                        tropomi_time = str(tropomi_dict['time'])[:5]
                        img_tropomi = folium.raster_layers.ImageOverlay(name=plume_name + '  |  ' + str(main_date.date()) + ' ' + tropomi_time + '  |  TROPOMI (S5p)',
                                                                        image='./FoliumMaps/Tropomi_Temp' + str(t) + '.png',
                                                                        bounds=[[cen_lat-delta_lonlat_tropomi, cen_lon-delta_lonlat_tropomi], [cen_lat+delta_lonlat_tropomi, cen_lon+delta_lonlat_tropomi]],
                                                                        origin='lower',
                                                                        show=False
                                                                       ).add_to(m)

                img = folium.raster_layers.ImageOverlay(name=plume_name + '  |  ' + str(main_date.date()) + ' ' + str(main_date.time())[0:5] + '  |  SENTINEL-3',
                                                        image='./FoliumMaps/S3_Temp.png',
                                                        bounds=[[cen_lat-delta_lonlat_S3, cen_lon-delta_lonlat_S3], [cen_lat+delta_lonlat_S3, cen_lon+delta_lonlat_S3]],
                                                        origin='lower',
                                                        show=False
                                                        ).add_to(m)


                if S2_dict: 
                    img = folium.raster_layers.ImageOverlay(name=plume_name + '  |  ' + str(main_date.date()) + ' ' + str(S2_dict['main_date'].time())[0:5] + '  |  SENTINEL-2',
                                                            image='./FoliumMaps/S2_Temp.png',
                                                            bounds=[[lonlat_bounds_S2[0][1], lonlat_bounds_S2[0][0]], [lonlat_bounds_S2[2][1], lonlat_bounds_S2[1][0]]],
                                                            colormap=cm.inferno,
                                                            origin='lower',
                                                            show=False
                                                           ).add_to(m)


            folium.Circle(radius=1000, location=[cen_lat, cen_lon], color='crimson', fill=False).add_to(m)
                    
            if tropomi_dicts:
                for t in range(len(tropomi_dicts)):
                    os.remove('./FoliumMaps/Tropomi_Temp' + str(t) + '.png')
            if S3_dict:
                os.remove('./FoliumMaps/S3_Temp.png')
            if S2_dict:
                os.remove('./FoliumMaps/S2_Temp.png')
    
    fig_folium_name += '.html'
    folium.LayerControl().add_to(m)
    m.save(fig_folium_name)
