from datetime import datetime as dt
import pandas as pd
import numpy as np

####################################
##         READ ME                ##
####################################

# The main dates can be specified in multiple ways through specific_main_dates
# 1): a list containing datetime date objects (excluding time stamps), [dt(y, m, d).date()]
specific_main_dates = [dt(2020, 10, 11).date(),
                       dt(2020, 10, 15).date()]
# 2): from a CSV file, with 'Days', 'Months', 'Years' as keys 
csv_file = '/data/maartenvn/Varon_Algeria/S2_detected_plume_dates.csv'
dates_df = pd.read_csv(csv_file)
days = dates_df['Days']
months = dates_df['Months']
years = dates_df['Years']
specific_main_dates = []
for i in range(len(days)):
    specific_main_dates.append(dt(years[i], months[i], days[i]))
# 3): 'All': all dates/observations within the all_data_dir directory
specific_main_dates = 'All'

####################################
source = '/data/maartenvn/Filtered_data/'
####################################

####################################
plume_name = 'Algeria'
####################################
cen_lon = 6.173
cen_lat = 31.864
delta_lonlat = 0.15
pixel_size = 0.001
all_data_dir = source + 'Algeria_Accident_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [#dt(2020, 1, 3).date(),
                       #dt(2020, 1, 4).date(),
                       dt(2020, 1, 5).date(),
#                        dt(2020, 1, 6).date(),
                       #dt(2020, 1, 7).date(),
                       #dt(2020, 1, 8).date(),
                       #dt(2020, 1, 9).date()
                        ]
algeria_accident_input = {'plume_name': plume_name,
                          'cen_lon': cen_lon,
                          'cen_lat': cen_lat,
                          'delta_lonlat': delta_lonlat,
                          'pixel_size': pixel_size,
                          'all_data_dir': all_data_dir,
                          'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Iraq'
####################################
cen_lon = 47.743
cen_lat = 30.274
delta_lonlat = 0.5
pixel_size = 0.001
all_data_dir = source + 'Iraq_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [#dt(2016, 7, 10).date(), 
                       #dt(2016, 7, 30).date(), 
                       #dt(2020, 8 , 8).date(),
                       #dt(2020, 10, 27).date(),
                       dt(2020, 7, 29).date(),
                       dt(2020, 9, 2).date()]
iraq_input = {'plume_name': plume_name,
              'cen_lon': cen_lon,
              'cen_lat': cen_lat,
              'delta_lonlat': delta_lonlat,
              'pixel_size': pixel_size,
              'all_data_dir': all_data_dir,
              'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Kazakhstan'
####################################
cen_lon = 59.0695
cen_lat = 49.3675
delta_lonlat = 0.15
pixel_size = 0.001
all_data_dir = source + 'Kazakhstan_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates =[dt(2021, 5, 14).date()]
kazakhstan_input = {'plume_name': plume_name,
                      'cen_lon': cen_lon,
                      'cen_lat': cen_lat,
                      'delta_lonlat': delta_lonlat,
                      'pixel_size': pixel_size,
                      'all_data_dir': all_data_dir,
                      'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Moscow_1'
####################################
cen_lon = 38.840
cen_lat = 56.685
delta_lonlat = 0.15
pixel_size = 0.001
all_data_dir = '/data/maartenvn/Moscow_Case/Moscow_1_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [dt(2021, 6, 18).date()]
moscow_1_input = {'plume_name': plume_name,
                  'cen_lon': cen_lon,
                  'cen_lat': cen_lat,
                  'delta_lonlat': delta_lonlat,
                  'pixel_size': pixel_size,
                  'all_data_dir': all_data_dir,
                  'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Moscow_2'
####################################
cen_lon = 38.489
cen_lat = 56.544
delta_lonlat = 0.15
pixel_size = 0.001
all_data_dir =  '/data/maartenvn/Moscow_Case/Moscow_2_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [dt(2021, 6, 18).date()]
moscow_2_input = {'plume_name': plume_name,
                  'cen_lon': cen_lon,
                  'cen_lat': cen_lat,
                  'delta_lonlat': delta_lonlat,
                  'pixel_size': pixel_size,
                  'all_data_dir': all_data_dir ,
                  'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Moscow'
####################################
# cen_lon = 38.66
# cen_lat = 56.61
delta_lonlat = 0.3
pixel_size = 0.001
# all_data_dir = source + 'Moscow_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'

cen_lon = 38.675
cen_lat = 56.675
all_data_dir = '/data/maartenvn/Moscow_Case/' + plume_name + '_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_newest/'

specific_main_dates = [dt(2021, 6, 18).date()] # [dt(y, m, d).date()]
moscow_both_input = {'plume_name': plume_name,
                      'cen_lon': cen_lon,
                      'cen_lat': cen_lat,
                      'delta_lonlat': delta_lonlat,
                      'pixel_size': pixel_size,
                      'all_data_dir': all_data_dir,
                      'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Permian'
####################################
cen_lon = -102.042349
cen_lat = 31.731678
delta_lonlat = 0.5
pixel_size = 0.001
all_data_dir = source + 'Permian_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [dt(2019, 3, 17).date()]
permian_input = {'plume_name': plume_name,
                  'cen_lon': cen_lon,
                  'cen_lat': cen_lat,
                  'delta_lonlat': delta_lonlat,
                  'pixel_size': pixel_size,
                  'all_data_dir': all_data_dir,
                  'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Sasol'
####################################
cen_lat = -26.553
cen_lon = 29.164
delta_lonlat = 0.15
pixel_size = 0.001
all_data_dir = '/data/maartenvn/Sasol_Case/Sasol_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'

csv_file = '/data/maartenvn/Sasol_Case/Plume_dates.csv'
dates_df = pd.read_csv(csv_file)
days = dates_df['Days']
months = dates_df['Months']
years = dates_df['Years']
specific_main_dates = []
for i in range(len(days)):
    specific_main_dates.append(dt(years[i], months[i], days[i]))

sasol_input = {'plume_name': plume_name,
                          'cen_lon': cen_lon,
                          'cen_lat': cen_lat,
                          'delta_lonlat': delta_lonlat,
                          'pixel_size': pixel_size,
                          'all_data_dir': all_data_dir,
                          'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Texas'
####################################
# cen_lon = -102.14; cen_lat = 31.40;  delta_lonlat = 0.5
cen_lon = -101.99; cen_lat = 31.45; delta_lonlat = 0.15
# cen_lon = -102.5; cen_lat = 31.05; delta_lonlat = 0.15
# cen_lon = -102.5; cen_lat = 31.95; delta_lonlat = 0.15
# cen_lon = -102.4; cen_lat = 31.85; delta_lonlat = 0.3
pixel_size = 0.001
all_data_dir = '/data/maartenvn/Texas/Texas_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [dt(2020, 9, 24).date(),
                       dt(2020, 10, 1).date(),
                       dt(2020, 10, 4).date(),
                       dt(2020, 10, 6).date(),
                       dt(2020, 10, 11).date(),
                       dt(2020, 10, 13).date(),
                       dt(2020, 10, 16).date(),
                       dt(2020, 11, 6).date(),
                       dt(2020, 11, 11).date(),
                       dt(2020, 12, 4).date(),
                       dt(2020, 12, 17).date()]
texas_input = {'plume_name': plume_name,
                          'cen_lon': cen_lon,
                          'cen_lat': cen_lat,
                          'delta_lonlat': delta_lonlat,
                          'pixel_size': pixel_size,
                          'all_data_dir': all_data_dir,
                          'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Dallas'
####################################
# cen_lon = -95.13; cen_lat = 32.0
cen_lon = -95.39; cen_lat = 32.07
delta_lonlat = 0.15
pixel_size = 0.001
all_data_dir = '/data/maartenvn/USA_Dallas/Dallas_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'
specific_main_dates = [dt(2021, 11, 29).date()]
dallas_input = {'plume_name': plume_name,
              'cen_lon': cen_lon,
              'cen_lat': cen_lat,
              'delta_lonlat': delta_lonlat,
              'pixel_size': pixel_size,
              'all_data_dir': '/data/maartenvn/' + all_data_dir,
              'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Varon_Algeria'
####################################
cen_lon = 5.9053
cen_lat = 31.6585
delta_lonlat = 0.5
pixel_size = 0.001
all_data_dir = source + 'Varon_Algeria_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'

option = 2
if option == 1:
    specific_main_dates = [dt(2020, 1, 29).date()]
elif option == 2:
    specific_main_dates = [dt(2019, 11, 20).date(),
#                            dt(2020, 1, 17).date(),
                           dt(2019, 11, 30).date(), 
                           dt(2019, 12, 5).date(), 
                           dt(2019, 12, 23).date(), 
                           dt(2019, 12, 25).date(), 
                           dt(2020, 1, 2).date(),
                           dt(2020, 1, 4).date(),
                           dt(2020, 1, 29).date(),
                           dt(2020, 2, 18).date(),
                           dt(2020, 2, 28).date(),
                           dt(2020, 3, 14).date()]
elif option == 3:
    csv_file = '/data/maartenvn/Varon_Algeria/S2_detected_plume_dates.csv'
    dates_df = pd.read_csv(csv_file)
    days = dates_df['Days']
    months = dates_df['Months']
    years = dates_df['Years']
    specific_main_dates = []
    for i in range(len(days)):
        specific_main_dates.append(dt(years[i], months[i], days[i]))
else:
    specific_main_dates = 'All'

varon_algeria_input = {'plume_name': plume_name,
                          'cen_lon': cen_lon,
                          'cen_lat': cen_lat,
                          'delta_lonlat': delta_lonlat,
                          'pixel_size': pixel_size,
                          'all_data_dir': all_data_dir,
                          'specific_main_dates': specific_main_dates}

####################################
plume_name = 'Varon_Turkmenistan'
####################################
cen_lon = 54.1977; cen_lat = 38.4939
# cen_lon = 53.9; cen_lat = 37.8
delta_lonlat = 0.5
pixel_size = 0.001
all_data_dir = source + 'Varon_Turkmenistan_' + str(np.around(cen_lon, decimals=4)) + '_' + str(np.around(cen_lat, decimals=4)) + '_interpolated_' + str(pixel_size) + '_' + str(delta_lonlat) + '/'

option = 1
if option == 1:
    specific_main_dates = [dt(2018, 7, 14).date(),
#                            dt(2019, 1, 15).date(),
#                            dt(2019, 2, 14).date(),
#                            dt(2019, 2, 24).date(),
#                            dt(2019, 3, 6).date(),
#                            dt(2019, 12, 21).date(),
#                            dt(2020, 3, 5).date(),
#                            dt(2020, 4, 29).date(),
#                            dt(2020, 10, 11).date()
                          ]
elif option == 2:
    csv_file = '/data/maartenvn/Varon_Turkmenistan/Dates_with_S3_plume.csv'
    dates_df = pd.read_csv(csv_file)
    days = dates_df['Days']
    months = dates_df['Months']
    years = dates_df['Years']
    specific_main_dates = []
    for i in range(len(days)):
        specific_main_dates.append(dt(years[i], months[i], days[i]))
elif option == 3:
    csv_file = '/data/maartenvn/Varon_Turkmenistan/PlumeDates.csv'
    dates_df = pd.read_csv(csv_file)
    days = dates_df['Days']
    months = dates_df['Months']
    years = dates_df['Years']
    specific_main_dates = []
    for i in range(len(days)):
        specific_main_dates.append(dt(years[i], months[i], days[i]))
else:
    specific_main_dates= 'All'

varon_turkmenistan_input = {'plume_name': plume_name,
                          'cen_lon': cen_lon,
                          'cen_lat': cen_lat,
                          'delta_lonlat': delta_lonlat,
                          'pixel_size': pixel_size,
                          'all_data_dir': all_data_dir,
                          'specific_main_dates': specific_main_dates}