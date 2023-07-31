#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import *
import pickle
import os,sys,glob
import csv
import numpy as np
from numpy.linalg import pinv
import matplotlib.pyplot as plt
#from typing import Callable
import pickle
import sys
from pylab import *
from datetime import datetime,date,timedelta


'''

Email: sudhanshu.pandey@jpl.nasa.gov
Institute: SRON Netherlands Institute for Space Research, Jet Propulsion Laboratory, Caltech, USA 


This script is part of the methane_detection code for methane plume detection from Multi spectral satellite  data. 
This code creates a lookup table for calculating the radiance of different atmospheric profiles. The radiance is calculated for different multi band satellites (Sentiel-2, Snetinel-3, Landsat)  and for different bands ("SWIR1", "SWIR2"). The script calculates the radiance for a range of zenith angles and different concentrations of methane in the atmosphere. The function `radianceCalc` calculates the radiance for a given satellite, band, solar zenith angle (sza), viewing zenith angle (vza), and methane concentration (del_omega). It uses absorption cross-sections for different gases (H2O, CO2, N2O, CO, CH4) and the number densities of these gases in different layers of the atmosphere. The radiance is then calculated by multiplying the solar irradiance at each wavelength by the transmission through the atmosphere, which is determined by the optical depth. The function `createLookUpTableFull` creates a lookup table for the radiance as a function of zenith angle and methane concentration. It loops over different satellites, zenith angles, and methane concentrations, and calculates the radiance for each combination. The results are stored in a dictionary and then saved to a pickle file.

'''

class ch4ret:
  def __init__(self, sza=0, vza= 0, method='MBMP',lat=31.6585,lng=5.9053,\
               area=4,satellite='S3',atmfile='atmosphere_midlatitudewinter.dat' ):
    """
    Intialize the object with main date, refdate, retrieval method, latitude and 
    longitude of location, study area (in km), windata satellite, TOA satellite and atmosphere profile
    """
    self.method=method # 'SBMP', 'MBSP' or 'MBMP'
    self.lat=lat
    self.lng=lng
    self.area=area # area in a km X a km
    self.atm_filename=atmfile
#    self.winddata_satellite=windsat # 'ERA5' or 'GFS' or 'GFS0P_Hourly'
    self.satellite=satellite # 'S2' or 'L8' or 'L7'
    self.delta=0.05
    self.coeff=np.zeros(6)
    self.markers=[[]]
    self.cdata={}
    self.sza= sza
    self.vza= vza



  def abCalc(self,band):
        """
        Function to generate .csv file containing absorption cross for each 
        molecule for each atmosphere layer in the selected profile.
        param band: Select band between B11 or B12
        
        """

        atm_filename=self.atm_filename

        if band=='SWIR1':
            hitran_filename='linelist2016_1600nm.dat'
        elif band=='SWIR2':
            hitran_filename='linelist2016_2100nm.dat'
        else:
            print('Invalid band selected, choose between B11 or B12')
          
        # Add the correct path to the HITRAN file
        data = np.loadtxt(hitran_filename)


        data_H2O=data[np.where(data[:,4]==1)]
        data_CO2=data[np.where(data[:,4]==2)]
        data_N2O=data[np.where(data[:,4]==4)]
        data_CO=data[np.where(data[:,4]==5)]
        data_CH4=data[np.where(data[:,4]==6)]

        tref=data[:,10][0] # reference temperature Tref [K]

        #HITRAN data for H2O molecule
        w0_H2O=data_H2O[:,0] # vacuum wavenumber [cm^-1]
        S0_H2O=data_H2O[:,1] # spectral line intensity at Tref [cm^-1/(mol cm^-2)]
        gamma_air_H2O=data_H2O[:,2] # air-broadened half-width half-maximum g at Tref [cm^-1 atm^-1]
        E_H2O=data_H2O[:,3] # lower-state energy [cm^-1] (unknown: 555.5555)
        mol_num_H2O=data_H2O[:,4] # molecule number
        mw_H2O=data_H2O[:,5] # molecular weight [u]
        beta_H2O=data_H2O[:,6] # coeff. for temp. dep. of ratio of total internal partition functions
        nair_H2O=data_H2O[:,7] # coeff. for temp. dep. of g
        delta_air_H2O=data_H2O[:,9] # air pressure-induced shift at Tref [cm^-1 atm^-1]

        
        #HITRAN data for CO2 molecule
        w0_CO2=data_CO2[:,0] # vacuum wavenumber [cm^-1]
        S0_CO2=data_CO2[:,1] # spectral line intensity at Tref [cm^-1/(mol cm^-2)]
        gamma_air_CO2=data_CO2[:,2] # air-broadened half-width half-maximum g at Tref [cm^-1 atm^-1]
        E_CO2=data_CO2[:,3] # lower-state energy [cm^-1] (unknown: 555.5555)
        mol_num_CO2=data_CO2[:,4] # molecule number
        mw_CO2=data_CO2[:,5] # molecular weight [u]
        beta_CO2=data_CO2[:,6] # coeff. for temp. dep. of ratio of total internal partition functions
        nair_CO2=data_CO2[:,7] # coeff. for temp. dep. of g
        delta_air_CO2=data_CO2[:,9] # air pressure-induced shift at Tref [cm^-1 atm^-1]

        #HITRAN data for N2O molecule
        w0_N2O=data_N2O[:,0] # vacuum wavenumber [cm^-1]
        S0_N2O=data_N2O[:,1] # spectral line intensity at Tref [cm^-1/(mol cm^-2)]
        gamma_air_N2O=data_N2O[:,2] # air-broadened half-width half-maximum g at Tref [cm^-1 atm^-1]
        E_N2O=data_N2O[:,3] # lower-state energy [cm^-1] (unknown: 555.5555)
        mol_num_N2O=data_N2O[:,4] # molecule number
        mw_N2O=data_N2O[:,5] # molecular weight [u]
        beta_N2O=data_N2O[:,6] # coeff. for temp. dep. of ratio of total internal partition functions
        nair_N2O=data_N2O[:,7] # coeff. for temp. dep. of g
        delta_air_N2O=data_N2O[:,9] # air pressure-induced shift at Tref [cm^-1 atm^-1]

        #HITRAN data for CO molecule
        w0_CO=data_CO[:,0] # vacuum wavenumber [cm^-1]
        S0_CO=data_CO[:,1] # spectral line intensity at Tref [cm^-1/(mol cm^-2)]
        gamma_air_CO=data_CO[:,2] # air-broadened half-width half-maximum g at Tref [cm^-1 atm^-1]
        E_CO=data_CO[:,3] # lower-state energy [cm^-1] (unknown: 555.5555)
        mol_num_CO=data_CO[:,4] # molecule number
        mw_CO=data_CO[:,5] # molecular weight [u]
        beta_CO=data_CO[:,6] # coeff. for temp. dep. of ratio of total internal partition functions
        nair_CO=data_CO[:,7] # coeff. for temp. dep. of g
        delta_air_CO=data_CO[:,9] # air pressure-induced shift at Tref [cm^-1 atm^-1]

        #HITRAN data for CH4 molecule
        w0_CH4=data_CH4[:,0] # vacuum wavenumber [cm^-1]
        S0_CH4=data_CH4[:,1] # spectral line intensity at Tref [cm^-1/(mol cm^-2)]
        gamma_air_CH4=data_CH4[:,2] # air-broadened half-width half-maximum g at Tref [cm^-1 atm^-1]
        E_CH4=data_CH4[:,3] # lower-state energy [cm^-1] (unknown: 555.5555)
        mol_num_CH4=data_CH4[:,4] # molecule number
        mw_CH4=data_CH4[:,5] # molecular weight [u]
        beta_CH4=data_CH4[:,6] # coeff. for temp. dep. of ratio of total internal partition functions
        nair_CH4=data_CH4[:,7] # coeff. for temp. dep. of g
        delta_air_CH4=data_CH4[:,9] # air pressure-induced shift at Tref [cm^-1 atm^-1]

        #Choose the correct atmosphere file with correct location, update atm_filename
        # in ch4ret.__init__ function  

        data_abs = np.loadtxt(atm_filename)
        lt=data_abs[:,0]
        p=data_abs[:,1] # pressure averaged over layer [hPa]
        p=p/1013.25  # pressure expressed in atm
        T=data_abs[:,2] # temperature averaged over molecules [K]

        # Atmosphere layer functions used for calculating tau_vertical
        def gammaL(tref,T,p,gamma_air,nair):
            gammaL=(tref/T)**(nair)*gamma_air*p
            return gammaL

        def gammaG(w0,T,mw):
            c3=3.581163241e-07   # in K^(−1/2)
            gammaG=c3*w0*(T/mw)**0.5
            return gammaG

        def rq(tref,T,beta):# rq:ratio of total internal partition functions
            rQ=(tref/T)**(beta)
            return rQ

        # temperature-corrected spectral line intensity
        def specLine(S0,E,w0,beta,tref,T):
            c2= 1.438776877 # in cm K
            S=S0*rq(tref,T,beta)*(exp(-c2*E*T**(-1))-exp(-c2*(E+w0)*T**(-1)))/\
                                  (exp(-c2*E/tref)-exp(-c2*(E+w0)/tref))
            return S

        # Voigt function (provided by Dr. Paul Tol)
        # Constants defined for Voigt function
        sqrtLn2 = np.sqrt(np.log(2.))
        recSqrtPi = 1./np.sqrt(np.pi)

        L24 = np.sqrt(24.)/2**0.25  # sqrt(N)/2^(1/4)
        a24 = np.array([ 0.000000000000e+00,
                        -1.5137461654527820e-10,  4.9048217339494231e-09,  1.3310461621784953e-09, -3.0082822751372376e-08,
                        -1.9122258508123359e-08,  1.8738343487053238e-07,  2.5682641345559087e-07, -1.0856475789744469e-06,
                        -3.0388931839363019e-06,  4.1394617248398666e-06,  3.0471066083229116e-05,  2.4331415462599627e-05,
                        -2.0748431511424456e-04, -7.8166429956141891e-04, -4.9364269012799368e-04,  6.2150063629501625e-03,
                        3.3723366855316406e-02,  1.0838723484566790e-01,  2.6549639598807695e-01,  5.3611395357291292e-01,
                        9.2570871385886777e-01,  1.3948196733791201e+00,  1.8562864992055403e+00,  2.1978589365315413e+00])
        w24polynom = np.poly1d(a24)

        def hum1wei24 (x,y):
	        """ Complex error function combining Humlicek's and Weideman's rational approximations:

	            |x|+y>15:  Humlicek (JQSRT, 1982) rational approximation for region I;
	            else:      J.A.C. Weideman (SIAM-NA 1994); equation (38.I) and table I.

	            F. Schreier, JQSRT 112, pp. 1010-1025, 2011:  doi: 10.1016/j.jqsrt.2010.12.010 """

	        # For safety only. Probably slows down everything. Comment it if you always have arrays (or use assertion?).
	        #if isinstance(x,(int,float)):  x = np.array([x])  # np.atleast_1d(x)

	        t = y - 1j*x
	        w = t * recSqrtPi / (0.5 + t*t)  # Humlicek (1982) approx 1 for s>15

	        if y<15.0:
		        mask = abs(x)+y<15.       # returns true for interior points
		        iz  = -t[np.where(mask)]  # returns small complex array covering only the interior region
		        # the following five lines are only evaluated for the interior grid points
		        lpiz = L24 + iz  # wL - y + x*complex(0.,1.)
		        lmiz = L24 - iz  # wL + y - x*complex(0.,1.)
		        recLmiZ  = 1.0 / lmiz
		        Z  = lpiz * recLmiZ
		        w24 = (recSqrtPi + 2.0*recLmiZ*w24polynom(Z)) * recLmiZ
		        # replace asympotic Humlicek approximation by Weideman rational approximation in interior center region
		        np.place(w, mask, w24)
	        return w

        def Voigt (vGrid, vLine=0.0, S=1.0, gammaL=1.0, gammaD=1.0):
            """ Voigt profile normalized to one, multiplied with line strength. """
            rgD = sqrtLn2 / gammaD
            x   = rgD * (vGrid-vLine)
            y   = rgD * gammaL
            vgt = recSqrtPi*rgD*S * hum1wei24(x,y).real  # Voigt profile
            return vgt
       
        # Calculate absorption cross section for each molecule and write it into csv files

        if self.satellite== "L8":
          if band=='SWIR1':
            lambda_l=1570; lambda_u=1650
          elif band=='SWIR2':
            lambda_l=2110; lambda_u=2290

        
        if self.satellite== "S2":
          if band=='SWIR1':
            lambda_l=1560 ;lambda_u=1660
          elif band=='SWIR2':
            lambda_l=2090  ;lambda_u=2290
          else:
            print('Invalid band selected, choose between B11 or B12')

        if self.satellite== "S3":
          if band=='SWIR1':
            lambda_l= 1613.40 -60.68/2 
            lambda_u=1613.40 +60.68/2 
          elif band=='SWIR2':
            lambda_l=  2255.70 - 50.15/2
            lambda_u= 2255.70 + 50.15/2
          else:
            print('Invalid band selected, choose between B11 or B12')

            

        wavenumber1=np.round((1e+07)/lambda_l,1)
        wavenumber2=np.round((1e+07)/lambda_u,1)
        delta=0.05
        w=np.arange(wavenumber2,wavenumber1,delta) # wavenumber grid
        wavelength=1e+07/w # wavelength grid 

        nLayer=p.shape[0]
        sigma_m=np.zeros([5, nLayer, w.shape[0]])

        # Absorption cross section calculation

        print('Calculating absorption cross section for H2O molecule:')
        lidx_first = np.where(w0_H2O > w[0] - 25)[0][0]
        lidx_last = np.where(w0_H2O < w[-1] + 25)[0][-1]
        for j in range(nLayer):
          #print("Calculating for atmospheric layer ",j)
          f_l = 0*w
          for k in range(lidx_first, lidx_last + 1):
            #print(k)
            gL = gammaL(tref,T[j],p[j],gamma_air_H2O[k],nair_H2O[k])
            gG = gammaG(w0_H2O[k],T[j],mw_H2O[k])
            S = specLine(S0_H2O[k],E_H2O[k],w0_H2O[k],beta_H2O[k],tref,T[j])
            vL = w0_H2O[k]+delta_air_H2O[k]*p[j]
            w_first = w0_H2O[k] - 25
            w_last = w0_H2O[k] + 25
        
            vG = w[(w > w_first) & (w < w_last)]
            f_abl =Voigt(vG,vLine=vL,S=S,gammaL=gL,gammaD=gG)
            f_l = f_l + np.concatenate([0*w[w <= w_first], f_abl , 0*w[w >= w_last]])
          sigma_m[0, j] = f_l
        sigma_H2O_pl=np.transpose(sigma_m[0,:,:])

        fil_name = '%s_absorption_cs_H2O_'%self.satellite+band
        data = sigma_H2O_pl
        data = data.tolist()
        with open(fil_name+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(data)

        print('Calculating absorption cross section for CO2 molecule:') 
        lidx_first = np.where(w0_CO2 > w[0] - 25)[0][0]
        lidx_last = np.where(w0_CO2 < w[-1] + 25)[0][-1]
        for j in range(nLayer):
          #print("Calculating for atmospheric layer ",j)
          f_l = 0*w
          for k in range(lidx_first, lidx_last + 1):
            #print(k)
            gL = gammaL(tref,T[j],p[j],gamma_air_CO2[k],nair_CO2[k])
            gG = gammaG(w0_CO2[k],T[j],mw_CO2[k])
            S = specLine(S0_CO2[k],E_CO2[k],w0_CO2[k],beta_CO2[k],tref,T[j])
            vL = w0_CO2[k]+delta_air_CO2[k]*p[j]
            w_first = w0_CO2[k] - 25
            w_last = w0_CO2[k] + 25
        
            vG = w[(w > w_first) & (w < w_last)]
            f_abl =Voigt(vG,vLine=vL,S=S,gammaL=gL,gammaD=gG)
            f_l = f_l + np.concatenate([0*w[w <= w_first], f_abl , 0*w[w >= w_last]])
          sigma_m[1, j] = f_l
        sigma_CO2_pl=np.transpose(sigma_m[1,:,:])
        fil_name = '%s_absorption_cs_CO2_'%self.satellite+band
        data = sigma_CO2_pl
        data = data.tolist()
        with open(fil_name+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(data)
        
        print('Calculating absorption cross section for N2O molecule:')
        lidx_first = np.where(w0_N2O > w[0] - 25)[0][0]
        lidx_last = np.where(w0_N2O < w[-1] + 25)[0][-1]
        for j in range(nLayer):
          #print("Calculating for atmospheric layer ",j)
          f_l = 0*w
          for k in range(lidx_first, lidx_last + 1):
            #print(k)
            gL = gammaL(tref,T[j],p[j],gamma_air_N2O[k],nair_N2O[k])
            gG = gammaG(w0_N2O[k],T[j],mw_N2O[k])
            S = specLine(S0_N2O[k],E_N2O[k],w0_N2O[k],beta_N2O[k],tref,T[j])
            vL = w0_N2O[k]+delta_air_N2O[k]*p[j]
            w_first = w0_N2O[k] - 25
            w_last = w0_N2O[k] + 25
        
            vG = w[(w > w_first) & (w < w_last)]
            f_abl =Voigt(vG,vLine=vL,S=S,gammaL=gL,gammaD=gG)
            f_l = f_l + np.concatenate([0*w[w <= w_first], f_abl , 0*w[w >= w_last]])
          sigma_m[2, j] = f_l
        sigma_N2O_pl=np.transpose(sigma_m[2,:,:])
        fil_name = '%s_absorption_cs_N2O_'%self.satellite+band
        data = sigma_N2O_pl
        data = data.tolist()
        with open(fil_name+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(data)

        print('Calculating absorption cross section for CO molecule:')
        lidx_first = np.where(w0_CO > w[0] - 25)[0][0]
        lidx_last = np.where(w0_CO < w[-1] + 25)[0][-1]
        for j in range(nLayer):
          #print("Calculating for atmospheric layer ",j)
          f_l = 0*w
          for k in range(lidx_first, lidx_last + 1):
            #print(k)
            gL = gammaL(tref,T[j],p[j],gamma_air_CO[k],nair_CO[k])
            gG = gammaG(w0_CO[k],T[j],mw_CO[k])
            S = specLine(S0_CO[k],E_CO[k],w0_CO[k],beta_CO[k],tref,T[j])
            vL = w0_CO[k]+delta_air_CO[k]*p[j]
            w_first = w0_CO[k] - 25
            w_last = w0_CO[k] + 25
        
            vG = w[(w > w_first) & (w < w_last)]
            f_abl =Voigt(vG,vLine=vL,S=S,gammaL=gL,gammaD=gG)
            f_l = f_l + np.concatenate([0*w[w <= w_first], f_abl , 0*w[w >= w_last]])
          sigma_m[3, j] = f_l
        sigma_CO_pl=np.transpose(sigma_m[3,:,:])

        fil_name = '%s_absorption_cs_CO_'%self.satellite+band
        data = sigma_CO_pl
        data = data.tolist()
        with open(fil_name+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(data)

        print('Calculating absorption cross section for CH4 molecule:') 
        lidx_first = np.where(w0_CH4 > w[0] - 25)[0][0]
        lidx_last = np.where(w0_CH4 < w[-1] + 25)[0][-1]
        for j in range(nLayer):
          #print("Calculating for atmospheric layer ",j)
          f_l = 0*w
          for k in range(lidx_first, lidx_last + 1):
            #print(k)
            gL = gammaL(tref,T[j],p[j],gamma_air_CH4[k],nair_CH4[k])
            gG = gammaG(w0_CH4[k],T[j],mw_CH4[k])
            S = specLine(S0_CH4[k],E_CH4[k],w0_CH4[k],beta_CH4[k],tref,T[j])
            vL = w0_CH4[k]+delta_air_CH4[k]*p[j]
            w_first = w0_CH4[k] - 25
            w_last = w0_CH4[k] + 25
        
            vG = w[(w > w_first) & (w < w_last)]
            f_abl =Voigt(vG,vLine=vL,S=S,gammaL=gL,gammaD=gG)
            f_l = f_l + np.concatenate([0*w[w <= w_first], f_abl , 0*w[w >= w_last]])
          sigma_m[4, j] = f_l
        sigma_CH4_pl=np.transpose(sigma_m[4,:,:])

        fil_name = '%s_absorption_cs_CH4_'%self.satellite+band
        data = sigma_CH4_pl
        data = data.tolist()
        with open(fil_name+'.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerows(data)

  def radianceCalc(self,del_omega,band, layer = 0):
      """
      Function to calculate spectral radiance integrated over given band 
      based on del_omega (mol/m2) added to the first layer of atmosphere

      param del_omega: methance column enhancement in mol/m2
      param band: string to select S2 band: 'B11' and 'B12'
      return spectral radiance integrated over chosen band range & optical depth
      of atmosphere molecules
      

      """

      atm_filename=self.atm_filename
      delta=self.delta

      if band=='SWIR1':
            E_filename='solar_irradiance_1600nm_highres_extended_sparse.dat'
      elif band=='SWIR2':
            E_filename='solar_irradiance_2100nm_highres_sparse.dat'
      else:
            print('Invalid band selected, choose between B11 or B12')
      
      sigma_H2O_pl= np.genfromtxt('%s_absorption_cs_H2O_'%self.satellite+band+'.csv', delimiter=',')
      sigma_CO2_pl= np.genfromtxt('%s_absorption_cs_CO2_'%self.satellite+band+'.csv', delimiter=',')
      sigma_N2O_pl= np.genfromtxt('%s_absorption_cs_N2O_'%self.satellite+band+'.csv', delimiter=',')
      sigma_CO_pl= np.genfromtxt('%s_absorption_cs_CO_'%self.satellite+band+'.csv', delimiter=',')
      sigma_CH4_pl= np.genfromtxt('%s_absorption_cs_CH4_'%self.satellite+band+'.csv', delimiter=',')
    
      data_abs = np.loadtxt(atm_filename)
      nH2O=data_abs[:,3] #column number density of H2O [cm^-2]
      nCO2=data_abs[:,4] #column number density of CO2 [cm^-2]
      nN2O=data_abs[:,6] #column number density of N2O [cm^-2]
      nCO=data_abs[:,7] #column number density of CO [cm^-2]\
      nCH4=data_abs[:,8] # column number density of CH4 [cm^-2]
      del_nCH4=del_omega*(6.023e+23)/10000
      nCH4[layer]=nCH4[layer]+del_nCH4

      nLayer=nH2O.shape[0]

      tau_vert = np.matmul(sigma_H2O_pl,nH2O[0:nLayer])+np.matmul(sigma_CO2_pl,nCO2[0:nLayer])\
      +np.matmul(sigma_N2O_pl,nN2O[0:nLayer])+np.matmul(sigma_CO_pl,nCO[0:nLayer])\
      +np.matmul(sigma_CH4_pl,nCH4[0:nLayer])

      optd_CO2=np.matmul(sigma_CO2_pl,nCO2[0:nLayer])
      optd_H2O=np.matmul(sigma_H2O_pl,nH2O[0:nLayer])
      optd_CH4=np.matmul(sigma_CH4_pl,nCH4[0:nLayer])

      if self.satellite== "L8":
          if band=='SWIR1':
            lambda_l=1570
            lambda_u=1650
          elif band=='SWIR2':
            lambda_l=2110
            lambda_u=2290
            
            
      if self.satellite== "S2":
          if band=='SWIR1':
            lambda_l=1560
            lambda_u=1660
          elif band=='SWIR2':
            lambda_l=2090
            lambda_u=2290
          else:
            print('Invalid band selected, choose between B11 or B12')

      if self.satellite== "S3":
          if band=='SWIR1':
            lambda_l= 1613.40 -60.68/2 
            lambda_u=1613.40 +60.68/2 
          elif band=='SWIR2':
            lambda_l=  2255.70 - 50.15/2
            lambda_u= 2255.70 + 50.15/2
          else:
            print('Invalid band selected, choose between B11 or B12')


      wavenumber1=np.round((1e+07)/lambda_l,1)
      wavenumber2=np.round((1e+07)/lambda_u,1)
      w=np.arange(wavenumber2,wavenumber1,delta)
      wavelength= 1e+07/w


      def f_young(za):
         za=radians(za)
         f=(1.002432*(cos(za))**2+0.148386*cos(za)+0.0096467)\
         /((cos(za))**3+0.149864*(cos(za))**2+0.0102963*cos(za)+0.000303978)
         return f
      
      A=0.1
      rair=f_young(self.sza)+f_young(self.vza)
      consTerm=A*cos(radians(self.sza))/pi
      #print("sza=",self.sza)
      tau_lambda=rair*tau_vert
      Edata = np.loadtxt(E_filename)
      
      E_lambda=np.zeros([w.shape[0]])
      w=np.round(w,2)
      
      for i in range(w.shape[0]):
          index=int(np.where(Edata[:,0]==w[i])[0])
          E_lambda[i]=Edata[index,1]

      
      L_lambda=consTerm*np.exp(-1*tau_lambda)*E_lambda
      self.L=L_lambda
      del_lambda=wavelength[0]-wavelength[1]
#      print (del_lambda, wavelength)
      T=np.sum(L_lambda)*del_lambda
      optd_data=np.transpose(np.vstack((wavelength,optd_H2O,optd_CO2,optd_CH4, tau_vert, E_lambda)))
      return T,optd_data
  def mCalc(self,del_omega,method):
      """
      Function to calcualte fractional change in spectral radiance based on 
      retrieval method selected

      param del_omega: methane column enhancements in mol/m2
      param method: Retrieval method- SBMP,MBSP,MBMP
      return fractional change in radiance, m
      
      # NOT NEEDED FOR GOOGLE COLAB
      """
      if method=='SBMP':
          T12_del,dump=self.radianceCalc(del_omega,'SWIR2')
          T12,dump=self.radianceCalc(0,'SWIR2')
          m=(T12_del-T12)/T12
      elif method=='MBSP' or 'MBMP':
          T12_del,dump=self.radianceCalc(del_omega,'SWIR2')
          T12,dump=self.radianceCalc(0,'SWIR2')
          T11_del,dump=self.radianceCalc(del_omega,'SWIR1')
          T11,dump=self.radianceCalc(0,'SWIR1')

          m=(T12_del-T12)/T12 - (T11_del-T11)/T11
      else:
          print('Select a valid retrieval method')
      return m
    
def f_young(za):
         za=radians(za)
         f=(1.002432*(cos(za))**2+0.148386*cos(za)+0.0096467)\
         /((cos(za))**3+0.149864*(cos(za))**2+0.0102963*cos(za)+0.000303978)
         return f
def giveamf(sza, vza):    
    rair=f_young(sza)+f_young(vza)
    return rair

def radianceCalc(del_omega,band, satellite, sza= 0, vza= 0, layer = 0, scaleb= 1.):
      delta=0.05
      atm_filename="atmosphere_midlatitudesummer.dat"
      if band=='SWIR1':
            E_filename='solar_irradiance_1600nm_highres_extended_sparse.dat'
      elif band=='SWIR2':
            E_filename='solar_irradiance_2100nm_highres_sparse.dat'
      
      sigma_H2O_pl= np.genfromtxt('%s_absorption_cs_H2O_'%satellite+band+'.csv', delimiter=',')
      sigma_CO2_pl= np.genfromtxt('%s_absorption_cs_CO2_'%satellite+band+'.csv', delimiter=',')
      sigma_N2O_pl= np.genfromtxt('%s_absorption_cs_N2O_'%satellite+band+'.csv', delimiter=',')
      sigma_CO_pl= np.genfromtxt('%s_absorption_cs_CO_'%satellite+band+'.csv', delimiter=',')
      sigma_CH4_pl= np.genfromtxt('%s_absorption_cs_CH4_'%satellite+band+'.csv', delimiter=',')
      data_abs = np.loadtxt(atm_filename)
      data_abs[0,3:9]= data_abs[0,3:9]*scaleb    
      nH2O=data_abs[:,3] #column number density of H2O [cm^-2]
      nCO2=data_abs[:,4] #column number density of CO2 [cm^-2]
      nN2O=data_abs[:,6] #column number density of N2O [cm^-2]
      nCO=data_abs[:,7] #column number density of CO [cm^-2]\
      nCH4=data_abs[:,8] # column number density of CH4 [cm^-2]
      del_nCH4=del_omega*(6.023e+23)/10000
      nCH4[layer]=nCH4[layer]+del_nCH4
      nLayer=nH2O.shape[0]

      tau_vert = np.matmul(sigma_H2O_pl,nH2O[0:nLayer] )+np.matmul(sigma_CO2_pl,nCO2[0:nLayer])\
      +np.matmul(sigma_N2O_pl,nN2O[0:nLayer])+np.matmul(sigma_CO_pl,nCO[0:nLayer])\
      +np.matmul(sigma_CH4_pl,nCH4[0:nLayer])

      optd_CO2=np.matmul(sigma_CO2_pl,nCO2[0:nLayer])
      optd_H2O=np.matmul(sigma_H2O_pl,nH2O[0:nLayer])
      optd_CH4=np.matmul(sigma_CH4_pl,nCH4[0:nLayer])

      if satellite== "S2":
          if band=='SWIR1':
            lambda_l=1560; lambda_u=1660
          elif band=='SWIR2':
            lambda_l=2090; lambda_u=2290
      if satellite== "S3":
          if band=='SWIR1':
            lambda_l= 1613.40 -60.68/2 ;            lambda_u=1613.40 +60.68/2 
          elif band=='SWIR2':
            lambda_l=  2255.70 - 50.15/2;            lambda_u= 2255.70 + 50.15/2
      wavenumber1=np.round((1e+07)/lambda_l,1)
      wavenumber2=np.round((1e+07)/lambda_u,1)
      w=np.arange(wavenumber2,wavenumber1,delta)
      wavelength= 1e+07/w

      
      A=0.1
      rair=f_young(sza)+f_young(vza)
      consTerm=A*cos(radians(sza))/pi
      tau_lambda=rair*tau_vert
      Edata = np.loadtxt(E_filename)

      E_lambda=np.zeros([w.shape[0]])
      w=np.round(w,2)
      
      for i in range(w.shape[0]):
          index=int(np.where(Edata[:,0]==w[i])[0])
          E_lambda[i]=Edata[index,1]

      
      L_lambda=consTerm*np.exp(-1*tau_lambda)*E_lambda
      L=L_lambda
      del_lambda=wavelength[0]-wavelength[1]
#      print (del_lambda, wavelength)
      T=np.sum(L_lambda)*del_lambda
      optd_data=np.transpose(np.vstack((wavelength,optd_H2O,optd_CO2,optd_CH4, tau_vert, E_lambda)))
      return T,optd_data


def mbmp(R12, R11, R12old, R11old): 
    """ Modified version of Multi band multi pass (Varon 2021) method to calculate radiance signal for methane.
    """
    R12= R12/median(R12)
    R11= R11/median(R11)
    R12old= R12old/median(R12old)
    R11old= R11old/median(R11old)
    delR1=np.divide((R12-R12old),R12old)   
    delR2=np.divide((R11-R11old),R11old)
    delR=delR1-delR2
    return delR





def giveMBMPgrid(omega,sat, sza, vza= 0):  
      omega_shape=  omega.shape
      R10, opd0 = radianceCalc(0,"SWIR1", sat,sza)
      R20, opd0 = radianceCalc(0,"SWIR2", sat,sza)
      omega_range= arange ( amin(omega),amax(omega), 0.25)
      mm= []
      for dell in omega_range:
        R2, opd0 = radianceCalc(dell,"SWIR2", sat,sza)
        R1, opd0 = radianceCalc(dell,"SWIR1", sat,sza)
        mfull= R2/R20 - R1/R10
        mm.append(mfull)
        print (dell )
      return omega_range, array(mm)  



def Omega2DelR(omega, satellite, sza , vza= 0, simple = True ):
    tamf = giveamf(sza,vza)
    shape_omega=    omega.shape 
    
    omega= omega.flatten()
    
    mdata= pickle.load(open('amf_mdata_poly_10_delr_to_omega.pkl', 'rb'))

    delr_data = mdata[satellite]["%2.1f"%tamf]
    delrr= []

    for ii, omm in enumerate(omega):
        ind = argmin(abs(mdata["omegas"] -  omm))
      #  print(delr_data[ind])
        delrr.append(delr_data[ind])
     #   print (omm,ii,  ind, delrr[ii], delr_data[ind],  mdata["omegas"][ind])
    delrr =  array(delrr).reshape(shape_omega)
    return delrr


def DelR2Omega(delr, satellite, sza,  method="MBMP" , vza= 0, simple = True ):
  """delr: single number or an 1d or 2d array.    
    Satellite: S2 or S3
    sza: solar zeinith angle 
    VZA: vewing zenith angle
    method: MBMP or MBSP or SBMP
    
    output
    omega: methane enhancement in mole/m^2
  """
  if simple:
    shape_omega=    delr.shape 
    delr= delr.flatten()
    mdata= pickle.load(open('amf_mdata_poly_10_delr_to_omega.pkl', 'rb'))
    tamf = giveamf(sza,vza)
    
#    print ("AMF", tamf, "%2.1f"%tamf)
    data = mdata[satellite]["%2.1f"%tamf]
    omega= zeros_like(delr)
    for ii, drr in enumerate(delr):
   
        ind = argmin(abs(data -  drr))
        omega[ii]=mdata["omegas"][ind]
      #  print (ii, drr, omega[ii])
    omega= omega.reshape(shape_omega)
    return omega
  else:
#        S2_full_mdata_poly_10_delr_to_omega.pkl
    mdata= pickle.load(open('S2_full_mdata_poly_10_delr_to_omega.pkl', 'rb'))
    tamf = (1/cos(radians (sza)) + 1/cos(radians (vza)  ))
    avamf= array(list(mdata[satellite][method].keys()))
    amf= avamf[argmin(abs(avamf- tamf))]
   # print ("AMF", tamf , amf)
    pol = poly1d(mdata[satellite][method][amf] )
    omega= pol(log(delr+1))
    return omega       




def giveOpticalDepth(satellite ):
    dat= {}#{'SWIR1': {} , "SWIR2": {} }
    for band in ["SWIR1", "SWIR2"]:
        obj= ch4ret(sza= 0, vza= 0, satellite = satellite) 
        R0, opd0 = obj.radianceCalc(0,band)
        R1, opd1 = obj.radianceCalc(1,band)
        E= opd0[:,-1]
        opt0=opd0[:,-2]
        dopt=opd1[:,-2] - opd0[:,-2]
        dat[band]= {"E" : E, "opt0" : opt0 , "dopt": dopt }
        
    return dat





def CreateAbsorptionCrossSections(sat= "S3"):
    os.chdir("./data_files")
    atmfile= ['atmosphere_europebackground.dat','atmosphere_subarcticsummer.dat','atmosphere_subarcticwinter.dat','atmosphere_tropical.dat', "atmosphere_midlatitudewinter.dat", "atmosphere_midlatitudesummer.dat","atmosphere_standard.dat" ]
    for aa in atmfile: 
        print (aa, os.path.exists(aa))
    method=['SBMP','MBMP']
  
    sza=np.arange(0,75,5)
    x1=np.arange(-2,5.1,0.1)
    x2=np.arange(5.5,15.5,0.5)
    obj= ch4ret(sza= 0, vza= 0, satellite = sat)
    x=np.concatenate((x1,x2))
    mdata={}
    for a in range(len(atmfile)):
      obj.atm_filename=atmfile[a]
      obj.abCalc('SWIR1')
      obj.abCalc('SWIR2')
    os.chdir("../")  

def plotOpticalDepths():
    """
    This function calculates optical depths for atmospheric gases and plots them using matplotlib. It's useful for visualizing the effect of gases on overall optical depth.
    """    
    os.chdir("./data_files")
    obj= ch4ret()
    obj.atm_filename='atmosphere_midlatitudewinter.dat'
    obj.abCalc('SWIR1')
    obj.abCalc('SWIR2')
    
    T11,optd_B11=obj.radianceCalc(0,'SWIR2', sza= 0, vza = 0)
    T11_,optd_B11_=obj.radianceCalc(1,'SWIR2', sza= 0, vza = 0)      
      
    plt.plot(optd_B11[:,0],optd_B11[:,1])
    plt.plot(optd_B11[:,0],optd_B11[:,2])
    plt.plot(optd_B11[:,0],optd_B11[:,3])
    plt.yscale('log')
    plt.xlabel('Wavelength')
    plt.ylabel('Optical depth (tau_vert)')
    plt.legend(['H2O','CO2','CH4'])
    plt.show()    
    os.chdir("../")  

def createLookupTables(sats= ["S3", "S2", "L8"]):
    """
    The function generates lookup tables for an atmospheric model based on different conditions. The main loop creates an instance of a class and calculates parameters for each combination of methods and air mass factors. The results are stored as a polynomial fit in a data dictionary and saved as a pickle file for later use.
    """

    
    os.chdir("./data_files")
    methods=['MBMP']
    szas=np.arange(0,75,10)
    x1=np.arange(-2,5.1,0.2)
    x2=np.arange(5.5,15.5,0.5)
    x=np.concatenate((x1,x2))

    start_time = datetime.now()
    data = {}

    total_methods = len(methods)
    total_amfs = len(np.arange(2,10,0.5))
    total_steps = total_methods * total_amfs

    current_step = 0

    for sat in sats:
        obj= ch4ret(sza= 0, vza= 0, satellite = sat)
        obj.atm_filename='atmosphere_midlatitudewinter.dat'
      
        data[sat]={}
        mdata= data[sat]
        
        for method in methods:
            mdata[method]={}
            for amf in np.arange(2,10,0.5):
                current_step += 1
                sza= 180*math.acos(1/ (amf-1))/np.pi
                mdata[method][amf]={}
                
                obj.sza= sza
                
                rads= []
                print (f"Processing: {sat}, {method}, sza: {sza:.2f}, amf: {amf}")
                print (f"Step {current_step} of {total_steps}")

                for xx in x: 
                    rads.append(obj.mCalc(xx,method))            

                rads= np.array(rads)      
                polyco= np.polyfit( np.log(rads+1),x, 10)
                mdata[method][amf]= polyco 

                print (f"Results for sza: {sza:.2f}, amf: {amf}, polynomial coefficients: {polyco}")
                print (f"Time elapsed: {datetime.now()- start_time}")

        pickle.dump(data, open(f'{sat}_full_mdata_poly_10_delr_to_omega.pkl', 'wb'))
        print(f"Data saved to {sat}_full_mdata_poly_10_delr_to_omega.pkl")
    os.chdir("../")  




if __name__ == '__main__':
    CreateAbsorptionCrossSections("L8")
    createLookupTables(["L8"]) # to create lookup tables for different atmospheric conditions
    plotOpticalDepths() # to make plot  of optical depths 
      
      
      
  
  
