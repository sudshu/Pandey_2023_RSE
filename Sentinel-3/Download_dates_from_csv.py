# Import modules
import os, shutil
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta as td
from suntime import Sun
from sentinelsat import SentinelAPI, geojson_to_wkt

####################################
##              Inputs            ##
####################################

head_folder = '/data/maartenvn/Texas/'

data_folder = head_folder + 'Texas/'

dates_file = head_folder + 'DownloadDates2.csv'

dates_to_download = pd.read_csv(dates_file)
days = dates_to_download['Days']
months = dates_to_download['Months']
years = dates_to_download['Years']

# Account information
username = 'maartensron'
password = 'maartensron'

longitude = -102.14
latitude = 31.40

# Data type
platformname = 'Sentinel-3'
# SL=SLSTR instrument, 1 = Level 1, RBT = Radiance and Brightness temperature
producttype = 'SL_1_RBT___'

####################################
##             DOWNLOAD           ##
####################################

# Login
api = SentinelAPI(username, password, 'https://scihub.copernicus.eu/dhus')
all_uuids_to_download = []

# Get the data information and id's to be downloaded
footprint = geojson_to_wkt({"type": "Feature","geometry": {"type": "Point","coordinates": [longitude, latitude]}})

for i in range(len(days)):
    day = days[i]
    month = months[i]
    year = years[i]
    
    date_start = dt(year, month, day)
    date_end = date_start + td(days=1)
    
    # Get the sunrise and sunset time in UTC for filtering the available datasets to only download sunlit observations
    sun = Sun(latitude, longitude)
    sunrise_utc = sun.get_sunrise_time(date_start).time()
    sunset_utc = sun.get_sunset_time(date_start).time()
        
    print('Period:  ', date_start, '  |  ', date_end)
    print('Location:  ', longitude, '  |  ', latitude)
    print('Sunrise:  ', sunrise_utc, '  |  Sunset:  ', sunset_utc)
    print()
        
    product = api.query(area=footprint,
                    date=(date_start, date_end),
                    platformname=platformname,
                    producttype=producttype
                       )
    product_df = api.to_dataframe(product)
    
    # Filter id's on overpass time, only consider daylight
    for uuid in product_df.index:
        info = api.get_product_odata(uuid)
        time = info['date'].time()

        if sunset_utc > sunrise_utc:
            if time > sunrise_utc and time < sunset_utc:
                all_uuids_to_download.append(uuid)
        else: 
            if time > sunrise_utc or time < sunset_utc:
                all_uuids_to_download.append(uuid)
    
print('Number of files to download: ', len(all_uuids_to_download))
print('Start downloading ... ')

# Download
api.download_all(all_uuids_to_download, directory_path = data_folder)

print(' ... Finished dowloading')

####################################
##               UNZIP            ##
####################################

print('Start unzipping ... ')

# Unzip the product files
for filename in os.listdir(data_folder):
    if filename.endswith('.zip'):
        shutil.unpack_archive(data_folder + filename, data_folder)
        os.remove(data_folder + filename)
        
print(' ... Finished unzipping')
