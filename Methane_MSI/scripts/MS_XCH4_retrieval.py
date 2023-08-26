

''' 
Sentinel-2, Landsat Methane retrieval code. This script contains the trunk of the code.
The code can be used as a module in a python notebook.

When importing this script as module in a python notebook, for autoreloading use:
     %load_ext autoreload
     %autoreload 2

Authors: Sudhanshu Pandey, Pratik Sutar, Paul Tol
Email: sudhanshu.pandey@jpl.nasa.gov
Institute: SRON Netherlands Institute for Space Research, Jet Propulsion Laboratory, Caltech, USA 

     
'''


import os,sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, Image
from pylab import * 
from datetime import datetime,date,timedelta
import ee, pickle
#ee.Authenticate() # only needed once
import collections; collections.Callable = collections.abc.Callable
ee.Initialize()



sat_info= {'S2':{"eepath":'COPERNICUS/S2','spacecraft':'SPACECRAFT_NAME','swir1':'B11',"swir2": "B12",\
                     "RGB" : ['B4', 'B3', 'B2'],'scaling':10000}, 
               "L8":{'eepath':"LANDSAT/LC08/C01/T1_TOA",'spacecraft':'SPACECRAFT_ID','swir1':'B6',"swir2":"B7",\
                 "RGB" : ['B4', 'B3', 'B2'],'scaling':1},
               "L7":{'eepath':"LANDSAT/LE07/C01/T1_TOA",'spacecraft':'SPACECRAFT_ID','swir1':'B5',"swir2":"B7",\
                     "RGB" : ['B3', 'B2', 'B1'],'scaling':1 },
               "L5":{'eepath':"LANDSAT/LE05/C01/T1_TOA",'spacecraft':'SPACECRAFT_ID','swir1':'B5',"swir2":"B7",\
                     "RGB" : ['B3', 'B2', 'B1'],'scaling':1 } }




class ch4ret:

  def __init__(self,main_date, ref_date,method='MBMP',lat=31.6585,lng=5.9053,\
               area=4,windsat='ERA5',TOAsat='S2',atmfile='atmosphere_midlatitudewinter.dat', verbose = False, case= 'temp', min_red_R= 0.5 ):
    """
    Intialize the object with maindate, refdate, retrieval method, latitude and 
    longitude of location, study area (in km), windata model, TOA satellite and atmosphere profile
    """
    self.main_date=main_date
    self.ref_date=ref_date # reference date
    self.method=method # 'SBMP', 'MBSP' or 'MBMP'
    self.lat=lat
    self.lng=lng
    self.area=area # area in a km X a km
    self.atm_filename=atmfile
    self.winddata_satellite=windsat # 'ERA5' or 'GFS' or 'GFS0P_Hourly'
    self.TOA_satellite=TOAsat # 'S2' or 'L8' or 'L7'
    self.delta=0.05
    self.coeff=np.zeros(6)
    self.markers=[[]]
    self.cdata={}
    self.thres=50 
    if self.TOA_satellite == 'S2': self.thres=10
    if self.TOA_satellite == 'L8': self.thres=50 
    if self.TOA_satellite == 'L8': self.thres=50 
    
    self.verbose = verbose
    self.area= area
    self.geom= ee.Geometry.Point((self.lng,self.lat))
    self.aoi= self.defineAOI(area)
    self.min_red_R= min_red_R      
    self.ref_index= ''
    self.main_index= ''
    self.case= case
    self.bad_indexs= []
    self.mdata= {}
    self.obs_ind= {}
    

  def defineAOI(self, area):
    "Defines Earth engine area of interest (AOI) object"
    buff=4900/10*area
    region = self.geom.buffer(buff).bounds().getInfo()['coordinates'][0]
    aoi = ee.Geometry.Polygon([[region[0],region[1],region[2],region[3]]], None, False)
    return aoi
    
  def check_obs_ind(self, cur_date):
  
    if len(self.obs_ind.keys())> 0:
      
      ddays = np.sort(np.array(list(self.obs_ind.keys()  )))
      if np.amax(ddays)< cur_date: return
      ii = where (ddays>= cur_date )[0][0]

      if  str(ddays[ii].date()) == str(self.main_date.date()): 
       
        return
      return self.obs_ind[ddays[ii]]
    return
  


  def getBandDataRange(self,start_date, end_date):
    """
    Function to search valid satellite scenes from "date" onwards

    """
    def saveData():
      filename= "output/" + self.case+'_banddata_raw'+ sat +  '.pkl'
      output = open(filename, 'wb')
      print ("# images",  count , dt.now()- st_time)
      
      pickle.dump(data, output)
      print ("saved to ", filename)
      output.close()
    from datetime import datetime as dt  
    aoi = self.aoi
    sdate = ee.Date(start_date)
    edate = ee.Date(end_date)    
    sat=self.TOA_satellite
    imgcol = ee.ImageCollection(sat_info[sat]["eepath"]).filterBounds(self.geom)\
    .filterDate(sdate, edate).sort('system:time_start')
    collist = imgcol.toList(imgcol.size().getInfo())
    data= {}
    st_time = dt.now()
    last_good_day = None
    cloudy_day= None
    for ii in range(collist.size().getInfo()):
        image = ee.Image(collist.get(ii))
        scene_date = datetime.utcfromtimestamp(image.getInfo()['properties']['system:time_start']/1000)
        day = str(scene_date.date())

        if day== last_good_day or cloudy_day ==day : continue
        band_arrs = image.sampleRectangle(region=aoi)
        index = image.getInfo()['properties']['system:index']
        
        if index in self.bad_indexs:          continue
        
        if self.verbose: print ("  " , str(scene_date.date())  , index, end= '')        
        try:
            red_data=np.array(band_arrs.get(sat_info[sat]["RGB"][0]).getInfo()) 
        except: 
          if self.verbose: print ("       Clipping issue")
          self.bad_indexs.append(index)
          continue


        if True:
              if sat=='S2':
                cdata=np.array(band_arrs.get('QA60').getInfo())  
                a=cdata>0
                CP = 100* (len(where(a.flatten())[0]) / len(cdata.flatten()))
                self.cloud_percentage= CP


              elif sat=='L8' or sat=='L7':
                self.cloud_percentage=    info = image.getInfo()['properties']['CLOUD_COVER']


                
#              if self.verbose: print( '     Cloud percentage : %2.2f '%(self.cloud_percentage))          

              if self.cloud_percentage > self.thres:
                if self.verbose : print ("       Cloudy")
                cloudy_day = day
                self.bad_indexs.append(index)          
                continue
             
        self.mdata[index]= {"reds" : np.hstack(red_data),  "date": scene_date  }
        
        self.obs_ind[scene_date]=  index
        
        if self.verbose: print (" .....  Sucessful observation")
        last_good_day= day 
        data[day]= {}
#        if self.verbose: print( "     ",   str(scene_date)[:19], ' Cloud percentage : %2.2f '%(self.cloud_percentage))
        R11,R12,rgb,acq_date,info=self.banddata(index = index)
        if sat=='S2':
        
           infos= { "sza" : info['MEAN_SOLAR_ZENITH_ANGLE'],
                    'index' : info['system:index'],                 'CF': info['CLOUDY_PIXEL_PERCENTAGE'] }
        elif sat=='L8' or sat=='L7':
            infos= { "sza" : 90- info['SUN_ELEVATION'],  'index' : info['system:index'],
                 'CF': info['CLOUD_COVER']                  }
           
        data[day]["R11"]= R11
        data[day]["R12"]= R12
        data[day]["rgb"]= rgb        
        data[day]['info']= infos
        count = len(data.keys())
#        print (ii,collist.size().getInfo())
        if count%10 ==0:
          saveData()

    saveData()

        
    return data 
    
    ndate= date+ timedelta(seconds=search_range*24*3600)
    print ("\n\n Warning : increasing search range as no valid observation found in days" , search_range)
    print (" New search date  %s  "%str(ndate.date()))

    
    if ndate> datetime.now():
      
      print ('No observation found till present')
      sys.exit()
      return 

    index = self.searchBandData(ndate)
    return index



  
  def searchBandData(self,date):
    """
    Function to search valid satellite scenes from "date" onwards

    """
    old_index = self.check_obs_ind(date)
    if not old_index  is None:
      return old_index
      
      
    aoi = self.aoi
    gdate = ee.Date(date)
  

    search_range= 60 #days forward
    sat=self.TOA_satellite
    imgcol = ee.ImageCollection(sat_info[sat]["eepath"]).filterBounds(self.geom)\
    .filterDate(gdate,gdate.advance(search_range,'day')).sort('system:time_start')
    collist = imgcol.toList(imgcol.size().getInfo())
    
    for i in range(collist.size().getInfo()):
        image = ee.Image(collist.get(i))
        scene_date = datetime.utcfromtimestamp(image.getInfo()['properties']['system:time_start']/1000)
        band_arrs = image.sampleRectangle(region=aoi)
        index = image.getInfo()['properties']['system:index']
        
        if index in self.bad_indexs:          continue
        
        if self.verbose: print ("  " , str(scene_date.date())  , index, end= '')        
        try:
            red_data=np.array(band_arrs.get(sat_info[sat]["RGB"][0]).getInfo()) 
        except: 
          if self.verbose: print ("       Clipping issue")
          self.bad_indexs.append(index)
          continue


        if True:
              if sat=='S2':
                cdata=np.array(band_arrs.get('QA60').getInfo())  
                a=cdata>0
                CP = 100* (len(where(a.flatten())[0]) / len(cdata.flatten()))
                self.cloud_percentage= CP


              elif sat=='L8' or sat=='L7':
                self.cloud_percentage=    info = image.getInfo()['properties']['CLOUD_COVER']


                
#              if self.verbose: print( '     Cloud percentage : %2.2f '%(self.cloud_percentage))          

              if self.cloud_percentage > self.thres:
                if self.verbose : print ("       Cloudy")
                self.bad_indexs.append(index)          
                continue
             
        self.mdata[index]= {"reds" : np.hstack(red_data),  "date": scene_date  }
        
        self.obs_ind[scene_date]=  index
        
        if self.verbose: print (" .....  Sucessful observation")

#        if self.verbose: print( "     ",   str(scene_date)[:19], ' Cloud percentage : %2.2f '%(self.cloud_percentage))          
        return index
    
    ndate= date+ timedelta(seconds=search_range*24*3600)
    print ("\n\n Warning : increasing search range as no valid observation found in days" , search_range)
    print (" New search date  %s  "%str(ndate.date()))

    
    if ndate> datetime.now():
      
      print ('No observation found till present')
      sys.exit()
      return 

    index = self.searchBandData(ndate)
    return index
 
  def banddata(self,date= None, index = None):
    """
    Function to extract TOA reflectance data from SWIR 1 and 2 band and store them
    ndarray. Also, extract RGB band data for maximum of 4x4 km area.

    param date: date for which band data is needed
    return np_arr_b11,np_arr_b12: SWIR 1 and SWIR 2 band TOA reflectance data
           rgb: Stacked matrix of RGB band data
           acq_date: Acquisition date for S2 image
           info: image data
    """  
    self.aoi = self.defineAOI(self.area)
    sat= self.TOA_satellite
   
    if index is None:
      index = self.searchBandData(date)
#      if index is None:
#        return 
   
    
    image = ee.Image(sat_info[sat]["eepath"]+'/'+index)  
    info = image.getInfo()['properties']
    band_arrs = image.sampleRectangle(region=self.aoi)
    spacecraft=info[sat_info[sat]["spacecraft"]]
    scene_date = datetime.utcfromtimestamp(info['system:time_start']/1000)

    acq_date=scene_date
    #print("Data acquisition date:%s from %s"%(scene_date,spacecraft))

    swir_data={}
    swir1_data=np.array(band_arrs.get(sat_info[sat]["swir1"]).getInfo())/sat_info[sat]["scaling"]
    swir2_data=np.array(band_arrs.get(sat_info[sat]["swir2"]).getInfo())/sat_info[sat]["scaling"]    
    
    rgb_band = band_arrs
    rgb_data1=np.array(rgb_band.get(sat_info[sat]["RGB"][0]).getInfo())/sat_info[sat]["scaling"]
    rgb_data2=np.array(rgb_band.get(sat_info[sat]["RGB"][1]).getInfo())/sat_info[sat]["scaling"]
    rgb_data3=np.array(rgb_band.get(sat_info[sat]["RGB"][2]).getInfo())/sat_info[sat]["scaling"]
    rgb_data=2*np.stack((rgb_data1,rgb_data2,rgb_data3),axis=2)
        
    return swir1_data,swir2_data,rgb_data,acq_date,info
  


  def make_cross(self, lat, lng, olat, olng):
    """
    Mark some location on the ouput images
    """
    import math
    n = len(olat)
    cx = []
    cy = []
    area=self.area

    r = 6371.0
    rs = r * math.cos(lat * math.pi/180)
    deg_lat = area/2/r * 180/math.pi
    deg_lng = area/2/rs * 180/math.pi

    min_lat = lat-deg_lat
    max_lat = lat+deg_lat
    min_lng = lng-deg_lng
    max_lng = lng+deg_lng

    for i in range(0,n):
        cx.append((olng[i]-min_lng)/(2*deg_lng))
        cy.append(1.0-(olat[i]-min_lat)/(2*deg_lat))
        
    return  cx,cy 


  
  def delRcalc(self, ref_date='2019-10-06',plot_option=True,calc_wind=True, verbose=True, olat=None, olng=None):
    """
    Function to calculate visualize methane plumes radiance signature deltaR based on retrieval
    method

    param ref_date: Reference date for methane retrieval 
    param print_option: figure boolean input: True or False
    param verbose: print boolean input: True or False
    """

    method=self.method

    #ref_date=self.ref_date
    main_date= self.main_date
    lat=self.lat
    lng=self.lng
    area=self.area
    sat=self.TOA_satellite
    if self.main_index== '':
      R11,R12,rgb,acq_date,info=self.banddata(main_date)
      self.main_index=info['system:index']
    else:
      R11,R12,rgb,acq_date,info=self.banddata(main_date, index = self.main_index )


    if self.ref_index== '':
     
      R11old,R12old,rgb_old,acq_date_old,info_old=self.banddata(ref_date )
      self.ref_index=info_old['system:index']
     
    else:
      R11old,R12old,rgb_old,acq_date_old,info_old=self.banddata(ref_date, self.ref_index)
 

    if str(acq_date.date() )== str(acq_date_old.date())  : 
        	acq_date_old = acq_date_old + timedelta(seconds = 24*3600)
        	self.ref_index = self.searchBandData (acq_date_old)
    R11old,R12old,rgb_old,acq_date_old,info_old=self.banddata(acq_date_old, self.ref_index)	
                

      
    self.main_date=acq_date;self.ref_date=acq_date_old

    if sat=='S2':
      self.sza1=info['MEAN_SOLAR_ZENITH_ANGLE']; self.sza2=info_old['MEAN_SOLAR_ZENITH_ANGLE']
      self.cloud_percentage=info['CLOUDY_PIXEL_PERCENTAGE']
      satellite1=info['SPACECRAFT_NAME'];satellite2=info_old['SPACECRAFT_NAME']
      self.cloud_percentage2=info_old['CLOUDY_PIXEL_PERCENTAGE']
      try:
        self.vza1=info['MEAN_INCIDENCE_ZENITH_ANGLE_B12']
        self.vza2=info_old['MEAN_INCIDENCE_ZENITH_ANGLE_B12']
      except:
        self.vza1=7.5;self.vza2=7.5  
    
      
    elif sat=='L8' or sat=='L7':
      self.sza1=90-info['SUN_ELEVATION']
      self.sza2=90-info_old['SUN_ELEVATION']
      self.cloud_percentage=info['CLOUD_COVER']
      self.cloud_percentage2=info_old['CLOUD_COVER']

      self.vza1=7.5;self.vza2=7.5
      satellite1=info['SPACECRAFT_ID'];satellite2=info_old['SPACECRAFT_ID']
        
#    self.main_date = datetime.utcfromtimestamp(info['system:time_start']/1000)
#    self.ref_date = datetime.utcfromtimestamp(info_old['system:time_start']/1000)
    
    if R11.shape>R11old.shape:
      m,n=R11old.shape[0],R11old.shape[1]
      R11=R11[:m,:n]
      R12=R12[:m,:n]
    elif R11.shape<R11old.shape:
      m,n=R11.shape[0],R11.shape[1]
      R11old=R11old[:m,:n]
      R12old=R12old[:m,:n]
    
    print('  Data acquisition date:%s from %s'%(self.main_date,satellite1))
    print('  Data acquisition date:%s from %s'%(self.ref_date,satellite2))


    
    if rgb.shape>rgb_old.shape:
      m,n=rgb_old.shape[0],rgb_old.shape[1]
      rgb=rgb[:m,:n,:]
    elif rgb.shape<rgb_old.shape:
      m,n=rgb.shape[0],rgb.shape[1]
      rgb_old=rgb_old[:m,:n,:]

    corr= np.corrcoef(np.hstack(R11old), np.hstack(R11))[0,1]
    corr_red= np.corrcoef(np.hstack(rgb[:,:,0]), np.hstack(rgb_old[:,:,0]))[0,1]
#    corr_g= np.corrcoef(np.hstack(rgb[:,:,1]), np.hstack(rgb_old[:,:,1]))[0,1]
#    corr_b= np.corrcoef(np.hstack(rgb[:,:,2]), np.hstack(rgb_old[:,:,2]))[0,1]
    self.R11= R11.copy()
    self.R12= R12.copy()
    self.R11old= R11old.copy()
    self.R12old= R12old.copy()     


    self.corr= corr
    self.STATUS='CLOUDLESS'
    CP=0
    
    
    if True:
      if method=='SBMP':
        c,y_int=self.lsReg(R12,R12old)
        delR=np.divide((c*R12-R12old),R12old)
        mainR=R12;refR=R12old
       
      elif method=='MBSP':
        c,y_int=self.lsReg(R12,R11)
        delR=np.divide((c*R12-R11),R11)
        mainR=R12;refR=R11

      elif method=='MBMP':
        c,y_int=self.lsReg(R12,R11)
        delR1=np.divide((c*R12-R11),R11)   
        c_old,y_int=self.lsReg(R12old,R11old)
        delR2=np.divide((c_old*R12old-R11old),R11old)
        delR=delR1-delR2
        self.delR1=delR1
        self.delR2=delR2
        mainR=R12;refR=R12old

      else:
        print('  Invalid method selected')  

    self.corr=corr
    self.corr_red=corr_red
    self.max_sigma=(np.amax(mainR)- np.median(mainR))/ np.std(mainR) 
    self.delR=delR

    if self.verbose:
#      print('Cloud percentage: %2.2f || STATUS: %s'%(CP,self.STATUS))
      print('  Current Date: SZA= %3.2f || VZA= %3.2f || eeImage cloud percentage= %3.2f'%(self.sza1,self.vza1,self.cloud_percentage))
      print('  Reference Date: SZA= %3.2f || VZA= %3.2f || eeImage cloud percentage= %3.2f'%(self.sza2,self.vza2,self.cloud_percentage2))
      print('  B11 Correlation= %1.2f || Max value in B12= %1.2f with %1.2f sigma || Red correlation= %1.2f'\
            %(corr,np.amax(mainR),(np.amax(mainR)- np.median(mainR))/ np.std(mainR),corr_red))
    self.delR=delR
    self.rgb= rgb
    self.rgb_old= rgb_old
    self.u10=0;self.v10=0; self.wind_speed=0;self.wind_dir=0
    self.refR= refR
    self.mainR= mainR
        
    if calc_wind: self.getWindSpeed()
    

    if plot_option:
        fig= self.plotDelRImage(corr, corr_red, olat, olng)
        self.fig_delR= fig
      
    
    

  def getWindSpeed(self):
    """
    wind speed at the observation time
    """
    try:
      u10,v10,wind_speed,wind_dir= self.windData(self.main_date)
      self.u10=u10;self.v10=v10;self.wind_speed=wind_speed;self.wind_dir=wind_dir
    except:
      print ("  Could not generate wind data from ERA5, switching to GFS" )
      try:
        self.winddata_satellite='GFS'
        u10,v10,wind_speed,wind_dir= self.windData(self.main_date)
        self.u10=u10;self.v10=v10; self.wind_speed=wind_speed;self.wind_dir=wind_dir
      except:
        print("  Could not generate wind data from both ERA5 and GFS0P")



    
  def plotDelRImage(self, corr, corr_red, olat= None, olng= None):
    fig= plt.figure(figsize=(20, 10))
    ax=plt.subplot(1,2,1)
    plt.sca(ax)
    im = ax.imshow(self.delR,cmap='inferno')
    ax.set_title('  Method: %s | Location: %6.3f \N{DEGREE SIGN}N ,%6.3f \N{DEGREE SIGN}E \n Main date: %s | Ref date: %s \
       \n WS= %4.1f m/s (%3.2f deg) | area=%dx%d  km$^2$ '%(self.method ,self.lat,self.lng, \
                                                str(self.main_date)[:16],str(self.ref_date)[:16], self.wind_speed, self.wind_dir,\
                                                self.area,self.area ),loc='center')
       
    ax.plot(self.delR.shape[0]/2 , self.delR.shape[1]/2, marker = 'x' , color = 'white')
    if olat!=None and olng!=None:
        n=len(olat)
        cx, cy= self.make_cross(self.lat, self.lng, olat, olng)
        for i in range(0,n):
            ax.plot(self.delR.shape[0]*cx[i] , self.delR.shape[1]*cy[i], marker = 'x' , color = 'white')
    im.set_clim((-0.15,0.15))
    #im.set_clim((0,0.2))
    q=ax.quiver(self.delR.shape[0]*0.9,self.delR.shape[1]*0.9,self.u10,self.v10,width=0.005,color='w')

    plt.gca().axes.get_yaxis().set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)

    axx=plt.subplot(2,4,8)
    im= axx.imshow(self.refR,cmap='gray' )
    #im.set_clim((0.3,0.8))
    axx.plot(self.delR.shape[0]/2 , self.delR.shape[1]/2, marker = 'x' , color = 'red')
    
    if olat!=None and olng!=None:
        n=len(olat)
        cx, cy= self.make_cross(self.lat, self.lng, olat, olng)
        for i in range(0,n):
            axx.plot(self.delR.shape[0]*cx[i] , self.delR.shape[1]*cy[i], marker = 'x' , color = 'red')
    
    axx.set_title('Control band (B11 if MBSP else B12old)')
    axx.axes.get_yaxis().set_visible(False)
    axx.axes.get_xaxis().set_visible(False)

    axx=plt.subplot(2,4,7)
    im= axx.imshow(self.mainR,cmap='gray' )
    #im.set_clim((0.3,0.8))
    axx.plot(self.delR.shape[0]/2 , self.delR.shape[1]/2, marker = 'x' , color = 'red')
    if olat!=None and olng!=None:
        n=len(olat)
        cx, cy= self.make_cross(self.lat, self.lng, olat, olng)
        for i in range(0,n):
            axx.plot(self.delR.shape[0]*cx[i] , self.delR.shape[1]*cy[i], marker = 'x' , color = 'red')
    axx.set_title('Main band B12 | B11 Correlation r: %6.2f'%corr )
    axx.axes.get_yaxis().set_visible(False)
    axx.axes.get_xaxis().set_visible(False)

    axx=plt.subplot(2,4,3)
    im= axx.imshow(self.rgb)
    #im.set_clim((0,0.5))
    axx.plot(self.rgb.shape[0]/2 , self.rgb.shape[1]/2, marker = 'x' , color = 'white')
    if olat!=None and olng!=None:
        n=len(olat)
        cx, cy= self.make_cross(self.lat, self.lng, olat, olng)
        for i in range(0,n):
            axx.plot(self.rgb.shape[0]*cx[i] , self.rgb.shape[1]*cy[i], marker = 'x' , color = 'white')
    axx.set_title('Main day RGB | Red Correlation r =%6.2f'%corr_red)
    axx.axes.get_yaxis().set_visible(False)
    axx.axes.get_xaxis().set_visible(False)

    axx=plt.subplot(2,4,4)
    im= axx.imshow(self.rgb_old)
    #im.set_clim((0,0.5))
    axx.plot(self.rgb.shape[0]/2 , self.rgb.shape[1]/2, marker = 'x' , color = 'white')
    if olat!=None and olng!=None:
        n=len(olat)
        cx, cy= self.make_cross(self.lat, self.lng, olat, olng)
        for i in range(0,n):
            axx.plot(self.rgb.shape[0]*cx[i] , self.rgb.shape[1]*cy[i], marker = 'x' , color = 'white')
    axx.set_title('Reference day RGB')
    axx.axes.get_yaxis().set_visible(False)
    axx.axes.get_xaxis().set_visible(False)
    plt.close(fig)
    return fig



  def ref_search(self,main_date, advance_range= 20 ,ref_date=None, ref_search_area=1, cthreshold= 0.85):
        """
        Function to return a  reference date for the main selected date based on maximum Red band
        correlation within a selected time period
        param advance_range: number of days before and after the main date for ref search 
              main_date: main date for which methane retrieval is done
              ref_date: Guess for ref_date (if correlation based search fails, ref_date will be returned as ref_date)
              threshold: Threshold for B11 correlation (default set at 0.6)
        return ref_date
        """


        print ("  **** Performing Reference day search **** ")
        dates=[];corrs=[]
        
#        self.aoi['refs'] = self.defineAOI(ref_search_area)
        sat=self.TOA_satellite
        
        if self.main_index == '':
          self.main_index = self.searchBandData(main_date)
          
        mdate= self.mdata[self.main_index]['date']
        self.main_date= mdate
        
        if ref_date is None:   ref_date = mdate+ timedelta(seconds = 22*3600)

          
        if self.ref_index == ''  :
          self.ref_index = self.searchBandData(ref_date)
          
#        print  (str(mdate.date()),  str(self.mdata[self.ref_index]['date'].date()) )
        if str(mdate.date())== str(self.mdata[self.ref_index]['date'].date()):
          
          ref_date = self.mdata[self.ref_index]['date'] + timedelta(seconds = 22*3600)
          self.ref_index = self.searchBandData (ref_date)
            
        red= self.mdata[self.main_index]['reds']
        red_old = self.mdata [self.ref_index]["reds"]
        red, red_old= checkShape(red, red_old)
        threshold=np.corrcoef(np.hstack(red), np.hstack(red_old))[0,1]

          
                        
        
        if self.verbose: print( "  Start ref day R threshold %3.2f "%threshold)
        if threshold > cthreshold:  return ref_date
        thress= [threshold]
        indexx= [self.ref_index]
        
        for ii, ind in enumerate (self.mdata.keys()):
          ddate= self.mdata[ind]['date'] 
          if ddate.date()==mdate.date(): continue
          if ddate + timedelta ( days= 60) <  mdate: continue
          red_old = self.mdata[ind]["reds"]
          red, red_old= checkShape(red, red_old)          
          threshold=np.corrcoef(np.hstack(red), np.hstack(red_old))[0,1]
          if self.verbose: print("   ", ii, ddate, " ref_day R %3.2f "%threshold)
          thress.append(threshold)
          indexx.append(ind)
          
        thress= array(thress)
        indexx= array(indexx)
        self.ref_index= indexx[argmax(thress)]
        ref_date = self.mdata[indexx[argmax(thress)]]['date']
        if thress.max() > cthreshold : 
          return ref_date
        else:
          cthreshold= thress.max()
        
        ee_main_date1=   ee.Date(str(self.main_date.date())).advance( 1, 'days')
        ee_main_date2=   ee.Date(str(self.main_date.date())).advance(advance_range+1, 'days')
        imgcol = ee.ImageCollection(sat_info[sat]["eepath"]).filterBounds(self.geom)\
                                                            .filterDate(ee_main_date1, ee_main_date2  ).sort('system:time_start')
        collist = imgcol.toList(imgcol.size().getInfo())

        if imgcol.size().getInfo() <  1 :
             if thress.max() < self.min_red_R:
                print (" best R of %2.2f is  too small. skipping "% ( corrs.max() )  )
                return
             else:
               return ref_date 
        

        
        print ("   Searching ahead of main_date ")
        
        date=[];corr=[]; indexs= []
#        cur_date=st_date
        scene_date = self.main_date + timedelta ( seconds = -5*3600 )
          
        for i in range(collist.size().getInfo()):
          image = ee.Image(collist.get(i))
          info = image.getInfo()['properties']
          new_scene_date= datetime.utcfromtimestamp(info['system:time_start']/1000)
          index = info['system:index']
          if new_scene_date.date() == self.main_date.date(): continue
          if self.main_index==index : continue
          
          if index in self.bad_indexs:          continue

          if self.verbose: print ("   " , str(new_scene_date.date())  , index, end= '')        
          
          if new_scene_date < scene_date +   timedelta(seconds= 2*3600):
            if self.verbose: print ("      skipping" )
            scene_date= new_scene_date
            continue
          else:
            scene_date= new_scene_date

          band_arrs = image.sampleRectangle(region=self.aoi)
          try:
             red_data=np.array(band_arrs.get(sat_info[sat]["RGB"][0]).getInfo()) 
            
#              test_temp=np.array(band_arrs.get(sat_info[sat]["swir1"][0]).getInfo()) 
          except:
              if self.verbose: print ("       Clipping issue")            

              continue
          if True:
               if sat=='S2':
                 cdata=np.array(band_arrs.get('QA60').getInfo())  
                 a=cdata>0
                 CP = 100* (len(where(a.flatten())[0]) / len(cdata.flatten()))
     #            if self.verbose: print( str(scene_date)[:19], ' Cloud percentage : %2.2f '%(CP))

               elif sat=='L8' or sat=='L7':
                 self.cloud_percentage=    info = image.getInfo()['properties']['CLOUD_COVER']

     
               if CP>self.thres:    
                   if self.verbose : print ("   Cloudy")
                   continue

          
          self.mdata[index]= {"reds" : np.hstack(red_data),  "date": scene_date  }

          redold=np.array(band_arrs.get(sat_info[sat]["RGB"][0]).getInfo())/sat_info[sat]["scaling"]
        
          cur_date=scene_date.date()
          
          # if red.shape>redold.shape:
          #       m,n=redold.shape[0],redold.shape[1]
          #       red=red[:m,:n]

          # elif red.shape<redold.shape:
          #       m,n=red.shape[0],red.shape[1]
          #       redold=redold[:m,:n]

              
          # To remove 'main date' as a choice for ref_date
          if cur_date != mdate.date():

            
            red, redold= checkShape(red, redold)
#            redold= check_translation(red, redold)
#            import pdb        ;    pdb.set_trace()
            cc= np.corrcoef(np.hstack(redold), np.hstack(red))[0,1]
            
            corrs+= [cc]
            dates+=[scene_date ]
            indexs+= [index]
            if self.verbose:print ("  corr %2.2f"%cc)
#          cur_date= datetime.fromisoformat(str(acq_date)) + timedelta(seconds = 22*3600)
#          cur_date=str(cur_date.date())
        if len(corrs)== 0:
          if self.verbose: print ('  !!!!  Ref_search range advance_range too small. Switching to default', ref_date)
          self.ref_index= ''
          return ref_date
          
        corrs=np.array(corrs)
        dates=np.array(dates)
        indexs= np.array(indexs)

        if corrs.max() > threshold:
          rdate=dates[np.where(corrs==corrs.max())][0]
          rindex=indexs[np.where(corrs==corrs.max())][0]


          if corrs.max() < self.min_red_R:
            print (" best R of %2.2f is  too small. skipping "% ( corrs.max() )  )
            return 
        
          if self.verbose: print ("\n  selected %s  corr %2.2f  "%(str(rdate)[:21], corrs.max() )  )
          
          self.ref_index= rindex  
          ref_date = rdate 
        print ("  **** Reference day search finished **** ")
        
        return ref_date


         
  def lsReg(self,dataB12,dataB11,plot_option=False):
    """
    Function for least square linear regression with slope output and 
    y-intercept set to 0

    param dataB12: main dataset
    param dataB11: reference dataset
    param plot_option: True or False
    return c,y_int: slope and y-intercept for the linear fit
    """
    R12=np.zeros([1000000])
    R11=np.zeros([1000000])
    for i in range(dataB11.shape[0]-1):
       R12=np.concatenate((R12,dataB12[i,:]))
       R11=np.concatenate((R11,dataB11[i,:]))

    c,y_int=np.polyfit(R12,R11,1)
    #print("slope=",c)
    #print("y-intercept=",y_int)
    R11_fitted=np.zeros(len(R11))
    for i in range(len(R12)):
        R11_fitted[i]=c*R12[i]+y_int
    if plot_option:
       plt.scatter(R12, R11, alpha=0.2)
       plt.plot (R12,R11_fitted,'r')
       plt.xlabel('$R_{12}$', fontsize=12)
       #plt.ylabel('$R_{11}$ or $R^{\'}_{12}$ ', fontsize=12)
       plt.ylabel('$R_{11}$ ', fontsize=12)
       plt.xlim((0, 1))
       plt.ylim((0, 1))
       plt.show()
    return c,y_int

  def windData(self,date):
    """
    Function to return u10 and v10 data from GFS: Global Forecast System
    u10:u_component_of_wind_10m_above_ground (m/sec)
    v10:v_component_of_wind_10m_above_ground (m/sec)

    param date: Date in 'yyyy-mm-dd hour-minute-second.milliseconds' format
    return u10 and v10 in m/sec, wind_speed and wind_dir (in m/sec and degrees resp.)
    """
    satellite=self.winddata_satellite
    if satellite == 'GFS':
      spath='NOAA/GFS0P25'
      u10_call='u_component_of_wind_10m_above_ground'
      v10_call='v_component_of_wind_10m_above_ground'
      advance_hour = 3.1
    elif satellite=='ERA5':
      
      spath='ECMWF/ERA5_LAND/HOURLY'
      u10_call='u_component_of_wind_10m'
      v10_call='v_component_of_wind_10m'
      advance_hour= 0.9
    elif satellite=='GFS0P_Hourly':
      u10,v10,wind_speed,wind_dir=geosFPatTimeLocation(self.lng, self.lat, date)
      return u10,v10,wind_speed,wind_dir
      
    
    date=ee.Date(date)
    geom = ee.Geometry.Point((self.lng,self.lat))
    region2 = geom.buffer(1000).bounds().getInfo()['coordinates']
     
    dataset = ee.ImageCollection(spath).filterBounds(geom).\
    filterDate(date.advance(-advance_hour,'hour'),date.advance(advance_hour,'hour')).sort('system:time_start').first()
    dataset=ee.Image(dataset)
    info = dataset.getInfo()['properties']
    index=info['system:index']
    index=spath+'/'+index

    date=datetime.utcfromtimestamp(info['system:time_start']/1000)
    
    img = ee.Image(index)
    region2=region2[0]

    aoi2 = ee.Geometry.Polygon([[region2[0],region2[1],region2[2],region2[3]]], None, False)
    data= img.sampleRectangle(region=aoi2)
    # Get individual band arrays.
    u10 =data.get(u10_call)
    u10=np.array(u10.getInfo())
    v10 = data.get(v10_call)
    v10=np.array(v10.getInfo())
    u10=u10[0][0];v10=v10[0][0];
    wind_speed=(u10**2+v10**2)**0.5
    wind_dir=degrees(atan2(v10,u10))
    if self.verbose:
      print ("  ", satellite , "found wind at %17.16s , %5.1f m/s (wind_direction= % 3.2f deg)"%(date, wind_speed,wind_dir ))
    return u10,v10,wind_speed,wind_dir


  def seriesDelr(self, start_date, end_date, ref_search_range= 30, continue_old= False, olat= None, olng= None ):
    """
    The function produces a time series of DelR images for the give period start_date and end_date
    """
    
    print ("!!!! Starting Series MS methane retrievals for the period %s--- %s   !!!! \n"%(str(start_date.date()), str(end_date.date())))
    import pickle
    advance_range=ref_search_range
    count=0;tcount=0; 
    filename= "output/" + self.case+'_plumedata'+'.pkl'
    print ("here")
    plumedata={}
    if continue_old:
    	if os.path.exists(filename):
             plumedata= pickle.load(open(filename, 'rb'))
             start_date = sort(list(plumedata))[-1]
             print ("switching the start date to %s"%start_date )
             start_date= datetime.strptime(start_date, "%Y-%m-%d")
             print (start_date)
    cur_date= start_date
    os.system ('mkdir -p %s/%s'%("output",self.case) )
    st_sr= datetime.now()
    
    while  cur_date <= end_date:
      date1 = str(cur_date.date())
      print ("\n>>>>>" ,date1,  "<<<<")
      self.main_index = self.searchBandData (cur_date)
      
      ref_date=self.ref_search(cur_date, ref_search_range)
      if ref_date is None:
        self.main_date= self.main_date + timedelta(seconds = 22*3600)
        cur_date= self.main_date
        continue
      
      self.delRcalc( ref_date,plot_option=True,  olat= None, olng= None)
      
      delR= self.delR;      fig= self.fig_delR
      self.main_index= ''
      if delR.all()!=0:
          del_omega = np.zeros_like(delR)
         # del_omega,fig3=obj.delr_to_del_omega()
         # del_omega,fig3=self.optimizer_v2()
        #  display(fig)
#          del_omega = np.zeros_like(delR)
          if self.wind_dir<0: self.wind_dir=360+self.wind_dir
          image_fname = 'output/%s/%s.png'%(self.case,self.TOA_satellite+'_delR_'+str(self.main_date)[:10])
          count+=1        
          fig.savefig(image_fname)
          print ("saving DelR image %s"%image_fname) 
          plumedata[str(self.main_date)[:10]]={'delR':delR,'del_omega':del_omega, "main_date": str(self.main_date) , 'ref_date':str(self.ref_date)[:10],  'wind_speed':self.wind_speed,\
                                              'wind_dir':self.wind_dir,'B11_corr':self.corr,'Red_corr':self.corr_red,'SZA':self.sza1,'VZA':self.vza1,\
                                              'SZA_ref':self.sza2,'VZA_ref':self.vza2  ,                                            
                                               "R11" : self.R11, "R12": self.R12 , "R12old": self.R12old, "R11old": self.R11old     }                                      

            
 
          if count%2 ==0:
            
            output = open(filename, 'wb')
            pickle.dump(plumedata, output)
            output.close()
            
          print("--------------------------------------------", end = '')
          print (" time from start: " , str(datetime.now()-st_sr)[:7])
#      if self.main_date
      self.main_date= self.main_date + timedelta(seconds = 22*3600);tcount+=1
      cur_date= self.main_date


    output = open( "output/"+ self.case+'_plumedata_'+self.TOA_satellite+'.pkl', 'wb')
    pickle.dump(plumedata, output)
    output.close()



  
  def foliumMap(self):
    """
    # main function for Folium 
    Function to return a interactive folium map for methane detection on a 100 x100 km2 area
    
    return
          delR_map: Folium map containing methane retrivals and TROPOMI detections
    """
    import folium    

    def add_ee_layer(self, eeImageObject, visParams, name):
         map_id_dict = ee.Image(eeImageObject).getMapId(visParams)
         folium.raster_layers.TileLayer(tiles = map_id_dict['tile_fetcher'].url_format,
                                        attr = "Map Data &copy; <a href='https://earthengine.google.com/'>Google Earth Engine</a>",
                                        name = name,
                                        overlay = True,control = True).add_to(self)
      # Add EE drawing method to folium.
    folium.Map.add_ee_layer = add_ee_layer


    main_date= self.main_date
    ref_date= self.ref_date

    
    sat=self.TOA_satellite
    if self.main_index== '':

      R11,R12,rgb,acq_date,info=self.banddata(main_date)
      self.main_index=info['system:index']
    else:
      R11,R12,rgb,acq_date,info=self.banddata(main_date, index = self.main_index )

    if self.ref_index== '':
#      self.ref_search()
      R11old,R12old,rgb_old,acq_date_old,info_old=self.banddata(ref_date)
      self.ref_index=info_old['system:index']
    else:
      R11old,R12old,rgb_old,acq_date_old,info_old=self.banddata(ref_date,  self.ref_index)


    img1=  ee.Image(sat_info[sat]["eepath"]+'/'+self.main_index)
    img2=  ee.Image(sat_info[sat]["eepath"]+'/'+self.ref_index)

    from IPython.core.display import display
    
    print("main_date: %s || ref_date  %s "%(str(acq_date)[:19], str(acq_date_old)[:19] ) )

    
    c,y_int=self.lsReg(R12,R11)
    c1,y_int1=self.lsReg(R12,R12old)
    c2,y_int2=self.lsReg(R12old,R11old)
    
    delR_Map = folium.Map(location=[self.lat, self.lng], zoom_start=12,tiles=None)
    vis = {'min': -0.15,'max': 0.15,'palette': ['FFFFB2', 'FED976', 'FEB24C','FD8D3C', 'FC4E2A', 'E31A1C', 'B10026']}
    #vis = {'min': -0.15,'max': 0.15,'palette':  ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718','74A901','66A000', '529400', '3E8601', '207401', '056201','004C00', '023B01', '012E01', '011D01', '011301']}
    folium.TileLayer(tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'\
                       ,attr= 'Tiles', name='ArcGIS').add_to(delR_Map)                    

    
    delR_MBSP = img1.expression('-1*(c*b12/b11 - 1)', {'b12': img1.select(sat_info[sat]["swir2"]),\
                                                       'b11': img1.select(sat_info[sat]["swir1"]),'c': c}) 
    delR_SBMP = img1.expression('-1*(c*b12/b11 - 1)', {'b12':img1.select(sat_info[sat]["swir2"]),\
                                                    'b11': img2.select(sat_info[sat]["swir2"]),'c': c1})
    delR_MBSP_old =img2.expression('-1*(c*b12/b11 - 1)', {'b12': img2.select(sat_info[sat]["swir2"]),\
                                                       'b11': img2.select(sat_info[sat]["swir1"]),'c': c2})
    delR_MBMP = img1.expression('del1-del2', {'del1': delR_MBSP,'del2': delR_MBSP_old})
    ratio = img1.expression('-1*((b12/b11)-(b12o/b11o))', {'b12': img1.select(sat_info[sat]["swir2"])\
                                                      ,'b11': img1.select(sat_info[sat]["swir1"]),\
                                                     'b12o': img2.select(sat_info[sat]["swir2"])\
                                                      ,'b11o': img2.select(sat_info[sat]["swir1"])})
    
    delR_Map.add_ee_layer(img1,{'bands': sat_info[sat]["RGB"], \
                              'min':0*sat_info[sat]["scaling"],'max':0.4*sat_info[sat]["scaling"]},'RGB')
    delR_Map.add_ee_layer(img1,{'bands': sat_info[sat]["swir2"], \
                              'min':0.1*sat_info[sat]["scaling"],'max':1*sat_info[sat]["scaling"]},'SWIR2')
    delR_Map.add_ee_layer(img2,{'bands': sat_info[sat]["RGB"], \
                              'min':0*sat_info[sat]["scaling"],'max':0.4*sat_info[sat]["scaling"]},'RGB_ref')
    delR_Map.add_ee_layer(img2,{'bands': sat_info[sat]["swir2"], \
                              'min':0.1*sat_info[sat]["scaling"],'max':1*sat_info[sat]["scaling"]},'SWIR2_ref')
    delR_Map.add_ee_layer(delR_SBMP, vis,'delR:SBMP')
    delR_Map.add_ee_layer(delR_MBSP, vis,'delR:MBSP')
    delR_Map.add_ee_layer(delR_MBMP, vis,'delR:MBMP')
    delR_Map.add_ee_layer(ratio, vis,'delR:Ratio')
    
    
    
    palette= ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
    date=ee.Date(str(acq_date)[:10])
    img  = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CH4').filter(ee.Filter.date(date, date.advance(1,"days"))).mean()
    tropomi_visParams = {'bands':['CH4_column_volume_mixing_ratio_dry_air'],'min':1700, 'max': 2000, "palette":palette }
    delR_Map.add_ee_layer(img, tropomi_visParams, "TROPOMI")
    
    if len(self.markers[0])==0:self.markers= Markers=[['Given Location'],[[self.lat,self.lng]]]
        
    for i in range(len(self.markers[0])):
      folium.Marker(location=self.markers[1][i],popup=self.markers[0][i]).add_to(delR_Map)
    
    from folium.plugins import MeasureControl , MousePosition
    delR_Map.add_child(folium.LayerControl())


    MousePosition().add_to(delR_Map)
    delR_Map.add_child(MeasureControl())
    #delR_Map = HTML('<iframe srcdoc="{}" style="float:left; width: {}px; height: {}px; display:inline-block; width: 100%; margin: 0 auto; border: 2px solid black"></iframe>'\
              # .format(delR_Map.get_root().render().replace('"', '&quot;'),1000,500))
    
    return delR_Map


  def plume_vis(self,del_omega,plot_option=False):
    """
    Function to visualize plume much clearly

    param del_omega: Methane column enhancement in mol/m2
    
    # NOT NEEDED FOR GOOGLE COLAB
    """
    test=del_omega
    test[np.where(test<0)]=0
    #index=np.where(test>=np.percentile(test,95))
    #data_index=np.transpose(np.vstack((index[0],index[1])))
    #mask=np.zeros(submatrix.shape)
    #for i in range(data_index.shape[0]):
      #mask[data_index[i][0],data_index[i][1]]=submatrix[data_index[i][0],data_index[i][1]]
    from scipy import ndimage

    test=ndimage.median_filter(test, size=(3,3))
    
    
    fig, ax = plt.subplots(figsize=(10,10))
    plt.imshow(test,cmap='RdBu_r')
    plt.title(' Method: %s | Location: %6.3f \N{DEGREE SIGN}N ,%6.3f \N{DEGREE SIGN}E \n Main date: %s | Ref date: %s \
       \n WS= %4.1f m/s (%3.2f deg) | area=%dx%d  km$^2$ '%(self.method ,self.lat,self.lng, \
                                                str(self.main_date)[:16],str(self.ref_date)[:16], self.wind_speed, self.wind_dir,\
                                                self.area,self.area ),loc='center')
    #cbar = plt.colorbar(shrink= 0.5)
    plt.plot(del_omega.shape[0]/2 , del_omega.shape[1]/2, marker = 'x' , color = 'black')
    #cbar.set_label('Methane column enhancement (mol/m2)', labelpad=10,rotation=270) 
    q=ax.quiver(delR.shape[0]*0.9,delR.shape[1]*0.9,self.u10,self.v10)
    plt.clim((-2,2))
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.close(fig)
    
    if plot_option:
        return test,fig
    else:
        return test,None

  def area_split(self,lat,lng,area,m):
    """
    Function to split study area into smaller chunks in order to access Earth Engine data

    param lng,lat: coordinates for the center
    param area: Study area size in km
    param m: Grid size input
    return regions: list of corner coordinates for each sub-region to define ee domain
           centers: list center cooridnates for each sub-region
    
    # NOT NEEDED FOR GOOGLE COLAB

    """
    from mpl_toolkits.basemap import Basemap
    self.lat=lat; self.lng=lng
  
    fig=plt.figure(figsize=(9, 3))
    map = Basemap(width=area*1000,height=area*1000,
                  resolution='h',projection='stere',
                  lat_0=lat,lon_0=lng)
  
    lons, lats, x, y = map.makegrid(m+1, m+1, returnxy=True)
    ax = fig.add_subplot(121)
    map.scatter(x, y, marker='o')
    map.drawcoastlines()
    plt.show()

    lat_array=lats[:,0]
    lon_array=lons[0]

    regions=[]; centers=[]
    for i in range(lat_array.shape[0]-1):
      for j in range(lon_array.shape[0]-1):
        coord=[[[lon_array[j],lat_array[i]],[lon_array[j],lat_array[i+1]],\
                [lon_array[j+1],lat_array[i+1]],[lon_array[j+1],lat_array[i]],\
                [lon_array[j],lat_array[i]]]]
        regions=regions+coord

    for i in range(len(regions)):
      centers=centers+[[0.5*(regions[i][0][0]+regions[i][2][0]),0.5*(regions[i][0][1]+regions[i][2][1])]]

    self.sub_regions=regions
    self.sub_regions_centers=centers
    self.sub_area=area/m

    return regions,centers

  def optimizer_v2(self,plot_option=False):

    """
    Function for retrieving time efficient methane column enhancements using lookup tables 
    for m data 
    Takes 3 seconds for solve

    return del_omega: Array containing methane enhancements in mol/m2

    """

    import pickle    
    pkl_file = open('lookup_tables/mdata_'+self.atm_filename.split('.',1)[0]+'.pkl', 'rb')
    mdata = pickle.load(pkl_file)
    pkl_file.close()
    delR=self.delR
    del_omega=np.zeros(delR.shape)
    if delR.all()==0: return np.zeros(delR.shape),None
    if self.sza1>70:self.sza1=70
    if self.sza2>70:self.sza2=70
    if self.method=="MBMP":
      delR=self.delR
      for i in range(delR.shape[0]):
        for j in range(delR.shape[1]):
          sza_key=5*round(self.sza1/5) # Rounding solar zenith angle to closest multiple of 5 degrees
          data=mdata[self.method][sza_key]
          del_omega[i,j]=min(data.keys(), key=(lambda k: (delR[i,j]-data[k])**2))
           
    else:
      for i in range(delR.shape[0]):
        for j in range(delR.shape[1]):
          sza_key=5*round(self.sza1/5) # Rounding solar zenith angle to closest multiple of 5 degrees
          data=mdata[self.method][sza_key]
          del_omega[i,j]=min(data.keys(), key=(lambda k: (delR[i,j]-data[k])**2))

   # del_omega=np.round(del_omega,2)

    fig, ax = plt.subplots(figsize=(10,10))
    plt.imshow(del_omega,cmap='RdBu_r')
    plt.title(' Method: %s | Location: %6.3f \N{DEGREE SIGN}N ,%6.3f \N{DEGREE SIGN}E \n Main date: %s | Ref date: %s \
       \n WS= %4.1f m/s (%3.2f deg) | area=%dx%d  km$^2$ '%(self.method ,self.lat,self.lng, \
                                                str(self.main_date)[:16],str(self.ref_date)[:16], self.wind_speed, self.wind_dir,\
                                                self.area,self.area ),loc='center')
    cbar = plt.colorbar(shrink= 0.5)
    plt.plot(delR.shape[0]/2 , delR.shape[1]/2, marker = 'x' , color = 'black')
    cbar.set_label('Methane column enhancement (mol/m2)', labelpad=10,rotation=270) 
    q=ax.quiver(delR.shape[0]*0.9,delR.shape[1]*0.9,self.u10,self.v10)
    plt.clim((-2,2))
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.close(fig)
    
    if plot_option:
        return del_omega,fig
    else:
      return del_omega,None





  def sourceRateCalc(self, plume_dir= None ,plume_width = 80,  plume_length=300, method= 'CSF', verbose = False, ): 	 
     if plume_dir  is None:   plume_dir = self.wind_dir  	 
     print (plume_dir)
     em , sig_em = sourceRateCalc (del_omega= self.del_omega, case= "", plume_dir= plume_dir, date=self.main_date, 
         plume_width = plume_width,  plume_length=plume_length, method= method, verbose = verbose)     
     return em, sig_em

		
# function to be used for emission quantification

# Plume Polar History (PPH) plot method

def printt(outm):
    sys.stdout.write("\r" + str(outm)) 
    sys.stdout.flush()

def polar_coord(delR):
    """
    Function to generate polar coordiantes from cartesian based on delR/del_omega matrix
    param delR: Input matrix either delR or del_omega
    return coord_r, coord_th polar coordiantes (r,theta)
    """
    center_pixel=[20*round(delR.shape[0]/2),-20*round(delR.shape[0]/2)]
    coord_x=np.zeros(delR.shape); coord_y=np.zeros(delR.shape)
    coord_r=np.zeros(delR.shape);coord_th=np.zeros(delR.shape)
    for i in range(int(delR.shape[0])):
        for j in range(int(delR.shape[1])):
            coord_x[i,j]=20*j-center_pixel[0]
            coord_y[i,j]=-20*i-center_pixel[1]
            coord_r[i,j]=sqrt(coord_x[i,j]**2+coord_y[i,j]**2)
            coord_th[i,j]=degrees(atan2(coord_y[i,j],coord_x[i,j]))
            
            if coord_th[i,j]<0: coord_th[i,j]=360+coord_th[i,j]
    return coord_r,coord_th



def plume_signal(delR,coord_r,coord_th):
    """
    Function to generate PPH plot
    param delR: Input matrix either delR or del_omega
          coord_r, coord_th: polar coordinates
    return signal,pdata: plume data averaged in terms of radial and angular coord
    """
    r=coord_r.flatten(); R=delR.flatten(); th=coord_th.flatten()
    from scipy import stats
    pdata=stats.binned_statistic_2d(r,th,R,bins=[8,18])
    signal=pdata[0][0:2].mean(axis=0)
    #signal=(signal-signal.mean())/signal.std()
    #plt.pcolormesh(pdata[2],pdata[1],pdata[0])
    plt.show()
    return signal,pdata
        
def PPHplotshow(fname = "Varon_2021_Algeria_plumedata.pkl"):
    pkl_file = open(fname, 'rb')
    plumedata = pickle.load(pkl_file)
    pkl_file.close()

    test=np.array([plumedata[i]['PPH_signal'] for i in plumedata])
    for i in range(test.shape[0]):
        test[i,:]=test[i,:]-np.median(test,axis=0)

    fig, ax = plt.subplots(figsize=(30,30))
    img=ax.imshow(test.T)
    yticks=list(np.linspace(0,17,18))
    ytick_labels=list(20*np.linspace(0,17,18))
    xticks=list(np.linspace(0,len(list(plumedata.keys()))-1,len(list(plumedata.keys()))))
    ax.set_xticks(xticks)
    ax.set_xticklabels(list(plumedata.keys()),rotation=90)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)
    ax.set_ylabel('Polar Coordinates')
    wind_dir=np.array([plumedata[i]['wind_dir'] for i in plumedata])/20
    plt.plot(wind_dir,lw=0,marker='x',color='r')
    #fig.colorbar(img)
    img.set_clim(vmin=0.5, vmax=1)
    plt.show()    

def checkGeosFPFile(cur_time):
        ''' Checks if geosFP file for aparticular time is there. If not download it and returns the file name'''

        #folder= '/scratch/shared/sud/meteo_data/geofp/'
        #os.chdir(folder)
        
        fname= cur_time.strftime( "GEOS.fp.asm.tavg1_2d_slv_Nx.%Y%m%d_%H30.V01.nc4")
        if not os.path.exists(fname):
               import wget
               print( "File not found...downloading from server")
               online_path = cur_time.strftime ('https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/das/Y%Y/M%m/D%d/')+ fname
               wget.download(online_path)
               
        return fname


def checkShape(red, red_old):
     red= red.flatten()
     red_old = red_old.flatten()
     print (red.shape, red_old.shape, end = " ")
     if not red.shape == red_old.shape:
          if len(red.shape) == 2:
            fshape1 = amin([red.shape[0], red_old.shape[0]]  )
            fshape2= amin([red.shape[1], red_old.shape[1]]  )
            red= red[:fshape1, :][:, :fshape2]
            red_old= red_old[:fshape1,:][:, :fshape2]
          if len(red.shape)== 1:
            fshape1 = amin([red.shape[0], red_old.shape[0]]  )
            red= red[:fshape1]
            red_old= red_old[:fshape1]
     print ("after ", red.shape, red_old.shape)
     return red, red_old

def geosFPatTimeLocation(  lon0, lat0, cur_time, verbose= False):
        import netCDF4 as nc4
        fname = checkGeosFPFile(cur_time)
        
        def getLonLonIndexes(lon,lat, lon0 , lat0):
          lat_ind = np.argmin(abs(lat-lat0))
          lon_ind = np.argmin(abs(lon-lon0))
          return lon_ind, lat_ind
        
        fid=nc4.Dataset(fname)
        lat= fid.variables['lat'][:]
        lon= fid.variables['lon'][:]
        lon_ind, lat_ind = getLonLonIndexes(lon,lat, lon0 , lat0)
        u10= fid.variables['U10M'][:][0,lat_ind , lon_ind]
        v10= fid.variables['V10M'][:][0, lat_ind, lon_ind]
        u50= fid.variables['U50M'][:][0,lat_ind , lon_ind]
        v50= fid.variables['V50M'][:][0, lat_ind, lon_ind]
        fid.close()

        wind_speed=(u10**2+v10**2)**0.5
        wind_dir=degrees(atan2(v10,u10))
        return u10, v10, wind_speed,wind_dir


# Source rate quantification
def makeSourceRectangle(x, y, dx, dy, wind_dir):
    # makes the source box over the methane enhancement/plume region 
    from shapely.geometry import Polygon
    plume = Polygon([ (x + dx, y), (x - dx, y), (x - dx, y + dy), (x + dx, y + dy) ])
                
    from shapely.affinity import rotate # check: https://shapely.readthedocs.io/en/latest/manual.html
    angl_corr = 90 if dx < dy else 0
    plume = rotate(plume, -(angl_corr - wind_dir), origin = [x,y])
    
    return plume

def makeTransects( data, lon0, lat0, wind_dir, alpha = 20, triangle_plume = False, no_of_transects = 15, plume_length=100, plume_width = 200,  verbose = False):
    # makes transects for the CSF method    
    from shapely.geometry import Point, LineString
    plume_coords = [lon0,lat0]
    alpha = radians(alpha) 
    P0 = Point(lon0, lat0)
    l = plume_length
    if verbose: print ( "plume length , plume_width" , l , plume_width)
    if True:
        # background transect
        wind_direction = radians(wind_dir)
        if verbose: print("\nWind dir value (degrees): {}.\n".format(wind_dir))
                
        dx = plume_width ; dy = l + 0.1 * l # where does it come from?
        x = P0.x; y = P0.y
        plume_box = makeSourceRectangle(x, y, dx, dy, wind_dir)
            
        transect_lines = []
        wsp_transects = []
        iy = 1
        for dl in linspace(0, l, no_of_transects): # so from 0.05 up to 0.7;
            A = Point(plume_coords[0] + dl*np.cos(wind_direction), plume_coords[1] + dl*np.sin(wind_direction))
            B = Point( A.x + l * np.cos( pi/2 + wind_direction ), A.y + l * np.sin( pi/2 + wind_direction ) )
            C = Point( A.x + l * np.cos( 3*pi/2 + wind_direction ), A.y + l * np.sin( 3*pi/2 + wind_direction ) )
            bL = LineString((B,C))

            if plume_box.intersects(bL):
                L = plume_box.intersection(bL) # shapely object;
                iy += 1
                transect_lines.append(L)
                        
    return transect_lines, plume_box

def makeIMEsourceRegion(lon0, lat0, wind_dir, alpha = 20, triangle_plume = False, extent_deg = 200, plume_width = 200, verbose = False):
     # makes a circle for the IME method 
    from shapely.geometry  import Polygon, Point 
    plume_coords = [lon0,lat0]
    alpha = radians(alpha) # degrees to radians;
    P0 = Point(lon0, lat0)
    l = extent_deg
            
    if True:
        wind_direction = radians(wind_dir) # degrees to radians;
        if verbose: print("\nWind dir value (degrees): {}.\n".format(wind_dir))
                   
        if triangle_plume: # this part is really not needed and never used, but was received with the initial syntax so it was kept in case it will be used in the future;
            P1 = Point(P0.x + l/cos(alpha)*np.cos(wind_direction-alpha), P0.y + l/cos(alpha)*np.sin(wind_direction-alpha))
            P2 = Point(P0.x + l/cos(alpha)*cos(wind_direction+alpha), P0.y + l/cos(alpha)*sin(wind_direction+alpha))
            P3 = Point(P0.x + l/cos(alpha)*cos(wind_direction+alpha), P0.y + l/cos(alpha)*sin(wind_direction+alpha))
            plume = Polygon([[P0.x, P0.y], [P1.x, P1.y], [P2.x, P2.y]])
                            
        else:
            # very similar to the makeTransects function from CSF method, same rectangle size;


            dx = plume_width ; dy = extent_deg
            print ("dx, dy", dx, dy)
            x = P0.x; y = P0.y
            plume = makeSourceRectangle(x, y, dx, dy, wind_dir) # used if I want to make IME region a square;
            #plume = P0.buffer(extent_deg) # used if I want to make IME region a circle;
            A_point_dir = Point(plume_coords[0] + (dy / 2)*np.cos(wind_direction), plume_coords[1] + (dy / 2)*np.sin(wind_direction))
                            
    return plume, A_point_dir

def makeBgBox(lon0, lat0, wind_dir = -10000, verbose = False):
    # the background values are taken from the box made here. It is is upwind of the source point  
    from shapely.geometry  import Polygon 
    dlat = dlon = 100 # give multiple of same value at once;
    wind_direction = radians(wind_dir)  # degrees to radians

    dll = abs(100) # size of what?: of the background box (half length of a side);
    dl_edge_point = 0.1 # distance between edge of box and source point of true emission;
    dl_bc = - (dl_edge_point + dll) # distance of center of box from source point of true emission;
            
    # an understanding of these next few lines would be needed?: done;
    bgclon = lon0 + dl_bc * np.cos(wind_direction)
    bgclat = lat0 + dl_bc * np.sin(wind_direction)
    latmax = bgclat + dll
    latmin = bgclat - dll
    lonmax = bgclon + dll
    lonmin = bgclon - dll

    if False:
        lonmax = bgclon + dll * 3
        lonmin = bgclon - dll * 3

    temp_poly = Polygon( [[lonmax, latmax], [lonmax, latmin] , [lonmin, latmin] , [lonmin, latmax]] )
    from shapely.affinity import rotate
    bg_box = rotate(temp_poly, -(90-wind_dir), origin = [bgclon, bgclat,])
#    print (temp_poly.exterior.xy)
    return bg_box

def giveErrorAreaMean(data_bg, areay, sigy):

    y = mean(data_bg * areay / mean(areay))
    error_y = (1 / sum( 1/sigy**2 ) )**0.5
    
    return y, error_y
        
def getPixcelsinPolygen_vd(poly_shape, pixels):
    
    vd = []; area = []
    
    for ii in arange(len(pixels)):
        pixel = pixels[ii] 
        if pixel.intersects(poly_shape):
            area.append(pixel.intersection(poly_shape).area)
            vd.append(ii)

    return array(vd), array(area)

def IMEMethod(pixels, data,  lon0, lat0, wind_dir, wind_speed, cur_date = None, thresh = 0.1, plume_width =200 , plume_length= 200, lim1_bg = 48, lim2_bg = 52, make_plot = False, verbose = False):
    extent_deg=   plume_length  
    background_box = makeBgBox(lon0, lat0, wind_dir, verbose = verbose)
    if verbose: print("\n\n starting IMEmethod")
    vd_bg, bg_area = getPixcelsinPolygen_vd(background_box, pixels)
    
    select_bg_values = data['xch4'][vd_bg][ (data['xch4'][vd_bg] >= percentile(data['xch4'][vd_bg], lim1_bg)) & (data['xch4'][vd_bg] <= percentile(data['xch4'][vd_bg], lim2_bg)) ]
    select_bg_area_values = bg_area[ (data['xch4'][vd_bg] >= percentile(data['xch4'][vd_bg], lim1_bg)) & (data['xch4'][vd_bg] <= percentile(data['xch4'][vd_bg], lim2_bg)) ]
    select_bg_error_values = data['dxch4'][vd_bg][ (data['xch4'][vd_bg] >= percentile(data['xch4'][vd_bg], lim1_bg)) & (data['xch4'][vd_bg] <= percentile(data['xch4'][vd_bg], lim2_bg)) ]
    bg_xch4, sig_bg_xch4 = giveErrorAreaMean(select_bg_values, select_bg_area_values, select_bg_error_values)
    
    plume, A_point_dir = makeIMEsourceRegion(lon0, lat0, wind_dir, triangle_plume = False, extent_deg = extent_deg, plume_width =plume_width,   verbose = verbose)
    # this function is similar with the function makeTransects from CSF method;
    vd_plume, source_area = getPixcelsinPolygen_vd(plume, pixels)

    # for IME method (for the windspeed function) I have to put st as input to the quantify function, unlike previously as each pixel is taken here, while the transects had a different function;
    em_est, sig_em_est, icount, ch4_estimates, wsp_list = quantifyIMEmission( data, pixels, wind_speed, bg_xch4, A_point_dir, plume, vd_plume, th_ppb_set = thresh, L_scale =None, make_plot = make_plot, verbose = verbose)

    rel_std_em = sig_em_est / em_est
    
    if make_plot:
        print ("IME emission estimate %5.2f +- %5.2f t/hr"%(em_est*3.6/1000 , sig_em_est *3.6/1000))
        plotMethod(data,  plume, background_box, lon0, lat0)
        title(" IME emission  %5.2f +- %5.2f t/hr \n %s"%(em_est*3.6/1000 , sig_em_est *3.6/1000, cur_date.date() ))
#        show()        
    return em_est* 3.6/1000, sig_em_est*3.6/1000 #, rel_std_em, icount, ch4_estimates, wsp_list

def quantifyIMEmission( data, tropomi_pixels, wind_speed, bg_xch4, A_point_dir, plume, vd_plume, th_ppb_set, L_scale= None,  make_plot = False, verbose = False):

    #    ''' calcualtes the IME emissions'''
 
    ch4_estimates = []
    enhancement = []
    wsp_list = []
    icount = []
    total_area = []
#    print ("sdsd", bg_xch4)

    xxsum= 0
    if True:
        sum_tau_area = total_plume_area = 0
        accepted_pixels = where( data['xch4'] >= th_ppb_set ) # each "data" dict key has a value for the number of pixels in the zone of interest, that we are plotting;
        print (  )
        for pixi in vd_plume: # iterate over the accepted pixels only -> save computation time;
            if pixi not in accepted_pixels[0]: continue
            # from lat, long degrees get real distance in km:
            area_km = tropomi_pixels[pixi].area
            
            total_plume_area += area_km # total area of pixels used so far;
            xxsum+=  data['xch4'][pixi]

            sum_tau_area += ( (data['xch4'][pixi] - bg_xch4)  * area_km )
            if make_plot: plot(tropomi_pixels[pixi].centroid.x, tropomi_pixels[pixi].centroid.y, zorder = 11, marker = 'x', color = 'black')

            icount.append(pixi)
        
        if L_scale is None:
            L_scale = sqrt(total_plume_area)
            print ("using plume area sqrt as L")
            
        else:
            print ("using circle radius as L")
        
        err= std(data['xch4'])/  (xxsum/len(icount)) 
#        print (std(data['xch4'] ) ,  mean(data['xch4'][vd_plume]) )
        ch4_estimates = (sum_tau_area * wind_speed / ( L_scale ))
        mean_emis = ch4_estimates
        error = np.nan
#        rest= set(arange(len(data["xch4"]))   ).difference(set(icount))
#        print (rest, rest)
#        if verbose: print("Result (g/sec): {:.2f} +/- {:.2f}".format(mean_emis, error))
#       print (list(rest))
#        import pdb; pdb.set_trace()        
        print (err)#, std(data['xch4'][rest]),  (xxsum/len(icount)) )

        return mean_emis,mean_emis* err, icount, ch4_estimates, wsp_list

def giveMidXY(del_omega, dll=20): # 20 m for S2
    del_omega= zeros_like(del_omega)
    mid_x= (del_omega+  dll*np.arange( del_omega.shape[1] )) +dll/2
    mid_y= (del_omega.transpose() +  dll*np.arange( del_omega.shape[0])) +dll/2 ;
    mid_y= np.transpose(mid_y[:,::-1])
    return mid_x, mid_y

def plotMethod( data,  plume, background_box, lon0, lat0):
    # makes maps for the IME and CSF methods 
#    plt.hexbin(data['x'], data['y'], data['xch4'],alpha = 0.5 ) ;  plt.colorbar(shrink = 0.5, label = 'del_omega (mole/m2)')
    plt.plot(lon0,lat0, marker = 'x')
    plot(*plume.exterior.xy, zorder = 10, color = 'purple', lw = 2, alpha = 1)
#    try:
    if not background_box is None :plot(*background_box.exterior.xy, zorder = 10, color = 'black', lw = 0.5, alpha = 1, label = 'Background methane area')
    plot(lon0, lat0, zorder = 11, marker = "x", color = 'white', markersize = 3, label = 'Source')
       #    xlabel('meters');        ylabel('meters')
      #    title('Plume analysis {} method (date: {}-{:02d}-{:02d}_{:02d}) \n Emission: {:.2f} $\pm$ {:.2f} (normalised) / Thresh: {:.1f}{} \n Wind Vel: {:.2f} $\pm$ {:.2f} m/s / Location: {}'.format(method, st.year, st.month, st.day, st.hour, em_est/10, sig_emi/10, thresh, units, wind_speed, chosen_source))#    Q10 = plt.quiver(data['lon'][::12], data['lat'][::12], data['uu0'][::12], data['vv0'][::12], zorder = 20, alpha = 0.7, scale = 200)#, color = 'yellow')
      # to make arrows white: color = 'yellow';
      #    gca().quiverkey(Q10, 0.0, 1.1, 5, r'5 m/s', labelpos = 'S', coordinates = 'axes')

def quantifyCSFluxes(data, tropomi_pixels, wind_speed, transect_lines, bg_xch4, sig_bg_xch4, plume, vd_plume, th_ppb_set, make_plot = False, verbose = False, skip_transects= 0):
    
    iy = len(transect_lines)
    ch4_estimates = []
    sig_ch4_estimates = []
    wsp_list = []
    icount = []
    # for tline in transect_lines[st_tr:]: # it start from st_tr (value is 0) up to the end (: was used);
    for tline in transect_lines[skip_transects:][ ::-1]:
                
        iy -= 1
        cme_dis = counter = lin_cover = 0
        sigma_temp = []
        threshold_ppb = []
    
        for pixi in vd_plume:
            
            if tropomi_pixels[pixi].intersects(tline):
                
                threshold_ppb.append(data['xch4'][pixi] - bg_xch4)
                # why multiplying by 100?: to put it in distance, [km];
                line_len = tropomi_pixels[pixi].intersection(tline).length 

                cme_dis = cme_dis + (data['xch4'][pixi] - bg_xch4) * line_len
                lin_cover += line_len
                # making sure more than 1 pixel intersec a transect;
                counter += 1

        if counter > 2:
            emi = cme_dis * wind_speed 
            ch4_estimates.append(emi)

            icount.append(iy+1)

                  
            if verbose: print( 'Transect {}, Emission (t/hr): {:.2f}, Pixels Intersected: {}'.format((iy+1), emi * 3600/1e6 ,  counter))
            if make_plot: # like for verbose, plot transects only if given the command;
                plot(*tline.xy, zorder = 10, color = 'red', linewidth = 1.5, label = 'Transects')

                if (iy+1) in [1, 5, 10, 15, 20]:
                       text(tline.xy[0][0] - 0.2, tline.xy[1][0], '%2i'%iy,  zorder = 10, fontsize = 5)
     
    # make sure there are no errors later with no ch4_estimates values or std that are 0 as there is only one value;
    mean_emis = mean(ch4_estimates) if len(ch4_estimates) > 0 else np.nan
    error = std(ch4_estimates) if len(ch4_estimates) > 1 else np.nan
            
#    if verbose: print("CSF Result (gm/sec): {:.2f} +/- {:.2f}".format(mean_emis, error))
    
    return mean_emis, error, icount, ch4_estimates, wsp_list

def CSFMethod( pixels, data, lon0, lat0, wind_dir, wind_speed, cur_date, thresh = 0, plume_width = 200,  plume_length=100,  no_of_transects= 10,  lim1_bg = 48, lim2_bg = 52, bg_xch4= None, make_plot = False, verbose = False , skip_transects=0):

    if bg_xch4 is None:
        background_box = makeBgBox(lon0, lat0, wind_dir, verbose = verbose)
        vd_bg, bg_area = getPixcelsinPolygen_vd(background_box, pixels)
    
        select_bg_values = data['xch4'][vd_bg][ (data['xch4'][vd_bg] >= percentile(data['xch4'][vd_bg], lim1_bg)) & (data['xch4'][vd_bg] <= percentile(data['xch4'][vd_bg], lim2_bg)) ]
        select_bg_area_values = bg_area[ (data['xch4'][vd_bg] >= percentile(data['xch4'][vd_bg], lim1_bg)) & (data['xch4'][vd_bg] <= percentile(data['xch4'][vd_bg], lim2_bg)) ]
        select_bg_error_values = data['dxch4'][vd_bg][ (data['xch4'][vd_bg] >= percentile(data['xch4'][vd_bg], lim1_bg)) & (data['xch4'][vd_bg] <= percentile(data['xch4'][vd_bg], lim2_bg)) ]
        bg_xch4, sig_bg_xch4 = giveErrorAreaMean(select_bg_values, select_bg_area_values, select_bg_error_values)
    else:
        bg_xhc4 = 0 ; sig_bg_xch4= 0.001
        background_box= None
    transects, plume=  makeTransects(data, lon0, lat0 , wind_dir, bg_xch4,False, no_of_transects,  plume_length , plume_width, verbose = verbose)
    # get the pixels that intersect and their area of intersection:
    vd_plume, source_area = getPixcelsinPolygen_vd(plume, pixels)

    em_est, sig_em_est, icount, ch4_estimates, wsp_list = quantifyCSFluxes(data, pixels, wind_speed, transects, bg_xch4, sig_bg_xch4, plume, vd_plume, th_ppb_set = thresh, make_plot = make_plot, verbose = verbose, skip_transects= skip_transects)
    
    rel_std_em = sig_em_est / em_est
    
    if make_plot:

        print ("CSF emission estimate %5.1f +- %5.1f t/hr"%(em_est*3.6/1000 , sig_em_est/1000))
        plotMethod(data,  plume, background_box, lon0, lat0)
        title("%s  %2.1f +- %2.1f t/hr"%(cur_date.date(), em_est*3.6/1000 , sig_em_est *3.6/1000,  ))
    return em_est *3.6/1000, sig_em_est *3.6/1000#, rel_std_em, icount, ch4_estimates

def sourceRateCalc(del_omega= None, case= "", date=datetime(1900,1,1) , 
                   plume_dir= 0, wind_speed = 2, location = None, sat= 'S2', 
                   plume_width = 80,  plume_length=300, method= 'CSF', verbose = False, 
                   ax = None, fname = None, make_plot= True , skip_transects= 0, dll = None):

    if del_omega is None:
        import pickle;
        dat= pickle.load(open(fname, 'rb'));
    if dll is None:    
      dll= 20 if sat == "S2" else 30 # satellite_pixel_size
      if sat == "S3" : dll= 100 #if  else 30 # satellite_pixel_size    

    
    del_omega = del_omega  * 16.04 # covert del_omega to grams/m2 methane SI units
    
    mid_x, mid_y= giveMidXY(del_omega, dll)
    
    from shapely.geometry import Polygon 
    data = {'x' : [], 'y':[] , 'xch4': [] , 'dxch4': []}
    pixel_list= []
    
    for ii in np.arange(del_omega.shape[0]):
            for jj in np.arange(del_omega.shape[1]):
                pixel= Polygon([(mid_x[ii,jj]-dll/2 , mid_y[ii,jj]-dll/2),   (mid_x[ii,jj]+dll/2 , mid_y[ii,jj]-dll/2) ,(mid_x[ii,jj]+dll/2 , mid_y[ii,jj]+dll/2),   (mid_x[ii,jj]-dll/2 , mid_y[ii,jj]+dll/2)])
                pixel_list.append(pixel)
                data['x'].append(mid_x[ii,jj])
                data['y'].append(mid_y[ii,jj])
                data['xch4'].append(del_omega[ii,jj])
                data['dxch4'].append(del_omega[ii,jj]/10.)

    for kk in data.keys(): data[kk]= array(data[kk])
#    threadhold = del_omega *16.04 #percentile(del_omega,95 )    
#    if verbose: print ("threadhold" , threadhold)
    if location is None:
        x0 , y0 = amax(data['x'])/2, amax(data['y'])/2
    else:
        x0 , y0 = location [0], location [1]
    if ax is None:
            figure(figsize= (6,5))
            plt.pcolormesh(mid_x, mid_y, del_omega*2.9/16.04, vmax = 2* std(del_omega*2.9/16.04)  ) # ; clim(1,30)
            
            ax= gca()
    else:
            sca(ax)
    
    if method.upper() == 'IME':em, sig_em= IMEMethod(pixel_list, data, x0,y0, wind_dir= plume_dir, wind_speed= wind_speed, cur_date= date , thresh= threadhold ,make_plot= make_plot ,  plume_width = plume_width, plume_length=plume_length,  verbose = True)

    if method.upper() == 'CSF': em, sig_em = CSFMethod( pixel_list, data, x0,y0, wind_dir= plume_dir , wind_speed= wind_speed,cur_date= date ,  plume_width = plume_width,  plume_length=plume_length , no_of_transects = 15, bg_xch4= 0, make_plot= make_plot ,  verbose = verbose, skip_transects= skip_transects)

    if make_plot:
        colorbar(shrink = 0.5, label = "XCH4 (ppm) %s-%s"%(sat, case))
        ylabel ("meters")
        xlabel("meters")
        show()
       
    return em , sig_em

def PointSeriesValues(  lng= 5.9053, lat= 31.6585, start_date= datetime(2019,1, 10), end_date= datetime (2020,1, 10), tag= "", sat="S2", verbose = False):
    ''' gives a pandas dataframe contaning all the dats with satellite overapasses''' 
    
    from ipygee import chart
    sat_info= {  'S2' : { "eepath" :'COPERNICUS/S2' , 'swir1' : 'B11', "swir2": "B12", "RGB" : ['B4', 'B3', 'B2']}, 
           "L7" : { 'eepath': "LANDSAT/LE07/C01/T1_SR" , 'swir1' : 'B5', "swir2": "B7" ,"RGB" : ['B3', 'B2', 'B1'] }}
    main_date= ee.Date(start_date)
    ref_date = ee.Date(end_date) 
    test_site = ee.Geometry.Point([ lng, lat])
    test_site= test_site.buffer(60)
    time_series= ee.ImageCollection(sat_info[sat]["eepath"]).filterBounds(test_site).filterDate(start_date, end_date)
    if sat == 'S2':
     chart_ts = chart.Image.series(**{
        'imageCollection': time_series, 
        'region': test_site,
        'scale': 10,
        'bands': ['B2', sat_info[sat]["swir1"], sat_info[sat]["swir2"],'B6'],
        'label_bands':[tag +"B2", tag + sat_info[sat]["swir1"], tag +sat_info[sat]["swir2"],tag+'B6'],
        'properties':['CLOUD_COVERAGE_ASSESSMENT'],
        'label_properties':[tag+ 'CLOUD_COVER'] })
    else:
         chart_ts = chart.Image.series(**{
        'imageCollection': time_series, 
        'region': test_site,
        'scale': 10,
        'bands': ['B2', sat_info[sat]["swir1"], sat_info[sat]["swir2"],'B6'],
        'label_bands':[tag +"B2", tag + sat_info[sat]["swir1"], tag +sat_info[sat]["swir2"],tag+'B6'],
        'properties':['CLOUD_COVER_LAND'],
        'label_properties':[tag+ 'CLOUD_COVER']  })
        
    if verbose : 
         print (chart_ts.dataframe)
        
    return chart_ts.dataframe
 
 
def add_ee_layer(m, eeImageObject, visParams, name=''):
    import folium
    map_id_dict = ee.Image(eeImageObject).getMapId(visParams)
    folium.raster_layers.TileLayer(tiles = ee.Image(eeImageObject).getMapId(visParams)['tile_fetcher'].url_format,
                                    attr = "Map",name = name, overlay = True,control = True).add_to(m)

def makeTropomiXCH4(lng, lat, date1 , date2= None, vrange= [1700, 2000]):
    from datetime import datetime
    from folium.plugins import MeasureControl
    if date2 is None:
    	date2= date1.advance(1,"days")
    m= folium.Map(location=[lat, lng], zoom_start=12, height=500,tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',attr= 'Tiles ')
    palette= ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
    date1=ee.Date(date1)
    img  = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CH4').filter(ee.Filter.date(date1, date2)).mean()
#     scene_date = datetime.utcfromtimestamp(img.getInfo()['properties']['system:time_start']/1000)
#     print (scene_date)
    tropomi_visParams = {'bands':['CH4_column_volume_mixing_ratio_dry_air'],'min':vrange[0], 'max': vrange[1], "palette":palette }
    add_ee_layer(m, img, tropomi_visParams, "Tropomi")
    folium.Marker(location=[lat, lng]).add_to(m)
    m.add_child(folium.LayerControl())
    # folium.LatLngPopup().add_to(m)
    m.add_child(MeasureControl())
    return m 
 
                                                                 
def simpleFoliumMap(lng= 5.9053, lat= 31.6585,  main_date=datetime(2019,10,20), end_date=datetime(2019,10,30) , sat= 'S2', zoom_start= 10 ):
    sat_info= {  'S2' : { "eepath" :'COPERNICUS/S2' , 'swir1' : 'B11', "swir2": "B12", "RGB" : ['B4', 'B3', 'B2']}, 
           "L7" : { 'eepath': "LANDSAT/LE07/C01/T1_SR" , 'swir1' : 'B5', "swir2": "B7" ,"RGB" : ['B3', 'B2', 'B1'] }}
    import folium
    from folium.plugins import MeasureControl
  #  print (main_date, end_date)
    main_date= ee.Date(main_date)
    end_date= ee.Date(end_date)
    if end_date is None: end_date= main_date.advance(10, "days") 

    test_site = ee.Geometry.Point([ lng, lat])
    
    m= folium.Map(location=[lat, lng], zoom_start=14, height=500,\
              tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',attr= 'Tiles ')
  
    
    def add_ee_layer(m, eeImageObject, visParams, name=''):
        
        map_id_dict = ee.Image(eeImageObject).getMapId(visParams)
        folium.raster_layers.TileLayer(tiles = ee.Image(eeImageObject).getMapId(visParams)['tile_fetcher'].url_format,
                                       attr = "Map",name = name, overlay = True,control = True).add_to(m)

    
    palette= ['black', 'blue', 'purple', 'cyan', 'green', 'yellow', 'red']
    img  =ee.ImageCollection(sat_info[sat]["eepath"]).filterDate(main_date, end_date).filterBounds(test_site).sort('system:time_start').first() ;
    scene_date = datetime.utcfromtimestamp(img.getInfo()['properties']['system:time_start']/1000)
    print (scene_date)
    name=  str(scene_date.date())[:10] + '_' +    sat+'_' +sat_info[sat]["swir2"]


    add_ee_layer(m, img, {'bands':[sat_info[sat]["swir2"]],'min':1000, 'max': 8000 },'SWIR2')
    add_ee_layer(m, img, {'bands':[sat_info[sat]["swir1"]],'min':1000, 'max': 8000 }, 'SWIR1')
    add_ee_layer(m , img, {'bands': sat_info[sat]["RGB"]  , 'min':1, 'max': 10000}       , 'RGB')

    
#    add_ee_layer(m, img.select([sat_info[sat]["swir2"]]).divide(img.select([sat_info[sat]["swir1"]])), {'min':0, 'max': 1} , "B12/B11")

    folium.Marker(location=[lat, lng]).add_to(m)

    m.add_child(folium.LayerControl())
    folium.LatLngPopup().add_to(m)
    m.add_child(MeasureControl())
    return m


def SeriesUsingPickleFile(case= None, pkl_file= None):
    ''' This function used the pickle file create in SeriesXCH4Imag and can remake DelR images with altered seting'''
    if case is None and pkl_file is None:
        print ('Not enough valid inputs. Please provide a valid case or a pickle file')
        return

    if pkl_file is None:        pkl_fil= 'output/%s_plumedata.pkl'%case


    from ms_locations import giveCaseParams
    case= 'kazakistan'
    lng, lat , start_date, ref_date, end_date = giveCaseParams(case)
      

    dat= pickle.load(open(pkl_file, 'rb'))
    for day_str in dat.keys():
      print ('processing %s'%day_str)
      main_date = datetime.strptime(day_str, "%Y-%m-%d")
      ref_date  =  datetime.strptime(dat[day_str]["ref_date"], "%Y-%m-%d")
      windsat='ERA5' # 'ERA5' or 'GFS' or 'GFS0P_Hourly'
      TOAsat='S2' # 'S2' or 'L8'  or 'L7'

      obj=ch4ret(main_date,ref_date,'MBMP',lat,lng,area,windsat,TOAsat)

      #ref_date=obj.ref_search(advance_range,date1,ref_date)
      obj.delRcalc(ref_date,verbose=True,plot_option=True, olat=olat, olng=olng)
      delR= obj.delR
      fig= obj.fig_delR 

    
      display(fig)
      image_fname = 'output/%s/%s.png'%(case,obj.TOA_satellite+'_delR_'+str(obj.main_date)[:10])
      fig.savefig(image_fname)
      
      return obj
      
      
      obj = dayXCH4Image(lng=lng, lat=lat, case= case, area = 1, main_date= main_date, ref_date=ref_date )
#      obj.
      return obj
#      fmap= obj.foliumMap(main_date, ref_date)
#      fmap.save
      image_fname = 'output/%s/%s.png'%(case,obj.TOA_satellite+'_delR_'+str(obj.main_date)[:10])

      obj.fig.savefig(image_fname)
      
      

def f_young(za):
         za=radians(za)
         f=(1.002432*(cos(za))**2+0.148386*cos(za)+0.0096467)\
         /((cos(za))**3+0.149864*(cos(za))**2+0.0102963*cos(za)+0.000303978)
         return f
def giveamf(sza, vza):    
    rair=f_young(sza)+f_young(vza)
    return rair

def fullMBMP2Omega(delr, satellite, sza, vza= 0 ):
    import pickle
    shape_omega=    delr.shape 
    delr= delr.flatten()
    mdata= pickle.load(open('data_files/test_srf_210104_delr_to_omega.pkl', 'rb'))
    tamf = giveamf(sza,vza)
 #   print ("AMF", tamf, "%2.1f"%tamf)
    ind_0 = where(mdata["omegas"]== 0 )
    dat = mdata[satellite]["%2.1f"%tamf]
    omega= zeros_like(delr)
    aa= dat["SWIR2"]/ dat["SWIR2"][0]
    bb= dat["SWIR1"]/ dat["SWIR1"][0]
    mbmp=  aa/bb-1
    for ii, drr in enumerate(delr):
        ind = argmin(abs(mbmp -  drr))
        omega[ii]=mdata["omegas"][ind]
    omega= omega.reshape(shape_omega)
    return omega

 

def dayXCH4Image(lng=5.9053, lat=31.6585, area=2, case= 'varon_2021_algeria', main_date=datetime(2019,10,20), ref_date=datetime(2019,10,6), make_plot= True, show_plot= True, verbose = False):
    ''' methane retrieval for a single day'''
    print ("!!!! Starting single day MS methane retrieval for day  %10s using ref  %10s  !!!!"%(main_date.date(), ref_date.date()) )
    olat=[lat];
    olng=[lng];

    advance_range=11
    windsat='ERA5' # 'ERA5' or 'GFS' or 'GFS0P_Hourly'
    TOAsat='S2' # 'S2' or 'L8'  or 'L7'

    obj=ch4ret(main_date,ref_date,'MBMP',lat,lng,area,windsat,TOAsat, verbose = verbose)

    ref_date=obj.ref_search(advance_range,main_date,ref_date,area/2. )
    
    obj.delRcalc(ref_date,plot_option=True)
    
    display(obj.fig_delR)
    return obj

  
    f_del_omega,fig3=obj.plume_vis(del_omega,plot_option=make_plot)
    obj.fig_xch4= fig2
    temp_fname = 'output/%s/%s.png'%(case,case+'_'+obj.TOA_satellite+'_delR_'+str(obj.main_date)[:10])
    fig.savefig(temp_fname)
    print ("...saving %s"%temp_fname)
    #fig3.savefig('S2/%s.png'%(case+'_'+obj.TOA_satellite+'_pm_'+str(obj.date1)[:10]))
#    sys.exit()    # srateCalc(obj.date1,del_omega,obj.wind_speed,obj.wind_dir,L=250)
    return obj



def calc_extent(lon,lat,dist):
    '''This function calculates extent of map
    Inputs:
        lat,lon: location in degrees
        dist: dist to edge from centre
    '''
    import cartopy.geodesic as cgeo
    dist_cnr = np.sqrt(2*dist**2)
    top_left = cgeo.Geodesic().direct(points=(lon,lat),azimuths=-45,distances=dist_cnr).base[:,0:2][0]
    bot_right = cgeo.Geodesic().direct(points=(lon,lat),azimuths=135,distances=dist_cnr).base[:,0:2][0]

    extent = [top_left[0], bot_right[0], bot_right[1], top_left[1]]

    return extent

def mbmpk(R12, R11, R12old, R11old): 
    R12= R12/nanmedian(R12)
    R11= R11/nanmedian(R11)
    R12old= R12old/nanmedian(R12old)
    R11old= R11old/nanmedian(R11old)
    delR1=np.divide(R12,R11)
  #  imshow(delR1); show()
    delR2=np.divide(R12old,R11old)
#    imshow(delR1); show()

#    delR1= check_translation(delR2,delR1 , makeplot= True, verbose = True)

#    imshow(delR1); show()
    delR=delR1-delR2
    print (" R", rrr(delR1, delR2))
    return delR


  

def rrr(aa,bb):
    # gives the correlation coefficient between two arrays

    aa= aa*bb/bb 
    bb= bb*aa/aa
    aa= aa.flatten()
    bb= bb.flatten()
    aa= aa[~isnan(aa)]
    bb= bb[~isnan(bb)]

    cor= corrcoef(aa, bb)[0,1]
    return round (cor, 3)  

def runSomething():
    case= 'varon_2021_algeria'
    os.system('mkdir -p output/%s'%case)
    lng=5.9053;    lat=31.6585
    start_date=datetime(2018,10,26)
    end_date=datetime(2020,5,10)    

    print ('running script')
    obj = dayXCH4Image(lng=lng, lat=lat, case= case, main_date= start_date, ref_date=start_date, verbose = False)  ;
    
    SeriesXCH4Image(lng=lng, lat=lat, case= case, start_date=start_date, ref_date= start_date, ed_date=end_date, verbose = True) ;
    print ("time taken: " , datetime.now()-dd) 
    
  

def giveWindRangeERA5( lng, lat, main_date, verbose = True):
    eu10s= []
    start_date=datetime(2018,6,29); ref_date= start_date
    obj=ch4ret(start_date,ref_date,'MBMP',lat,lng,area= 2,case= "",TOAsat= 'S2', verbose = True, min_red_R= 0.7 )
    for hour in arange(-2,3,1):
        u10,v10,wind_speed,wind_dir= obj.windData(main_date+ timedelta(seconds= 3600.*hour) )
        if verbose: print (wind_speed,wind_dir)
        eu10s.append(wind_speed)
    wind_e = std(eu10s)
    return mean(eu10s) , wind_e


if __name__ == "__main__":
    runSomething()
    sys.exit()

    
    
    
