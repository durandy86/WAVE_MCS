#!/usr/bin/env python
# coding: utf-8

##########################################################################
# Histogrammes permettant de choisir les variables les plus pertinentes  #
# pour ensuite les utiliser dans les methodes de clustering              #
##########################################################################

### Importation des bibliotheques
import sys
import matplotlib as mpl
import time
from math import *
import numpy as np 
import gzip
import subprocess

from struct import unpack
from struct import *

import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.gridspec as gridspec
import gzip

from datetime import *
from math import *
from mpl_toolkits.mplot3d import Axes3D
from xml.dom import minidom    
from random import randint
# from datetime import datetime

##################################################################################
#
#  Class convective system_IntParameters   : 
#
#
##################################################################################

class MCS_IntParameters(object): 
    def __init__(self):
        self.label = 0
        self.qc_MCS = 0
        self.duration = 0
        self.nbslot = 0
        self.Utime_Init = 0
        self.lonInit = 0
        self.latInit = 0
        self.localtime_Init = 0
        self.Utime_End = 0
        self.lonEnd = 0
        self.latEnd = 0
        self.localtime_End = 0
        self.vavg = 0
        self.dist = 0
        self.lonmin = 0
        self.lonmax = 0
        self.latmin = 0
        self.latmax = 0
        self.tbmin = 0
        self.surfmaxPix_235K = 0
        self.surfmaxkm2_220K = 0
        self.surfmaxkm2_210K = 0
        self.surfmaxkm2_200K = 0
        self.surfcumkm2      = 0


class MCS_Lifecycle(object):

    def __init__(self):
        self.qc_im				= []
        self.tbmin				= []
        self.tbavg_235				= []
        self.tbavg_208				= []
        self.tbavg_200				= []
        self.tbavg_90th				= []
        self.surfPix_235K		= []
        self.surfPix_210K		= []
        self.surfKm2			= []
        self.Utime				= []
        self.Localtime			= []
        self.lon				= []
        self.lat				= []
        self.x					= []
        self.y					= []
        self.velocity			= []
        self.semiminor_220K		= []
        self.semimajor_220K		= []
        self.orientation_220K	= []
        self.excentricity_220K	= []
        self.semiminor_235K		= []
        self.semimajor_235K		= []
        self.orientation_235K	= []
        self.excentricity_235K	= []
        self.surfkm2_235K		= []
        self.surfkm2_220K		= []
        self.surfkm2_210K		= []
        self.surfkm2_200K		= []     

def load_TOOCAN(FileTOOCAN):

    lunit = open(FileTOOCAN,'r')

    #
    # Read the Header
    ##########################
    header1 = lunit.readline()
    header2 = lunit.readline()
    header3 = lunit.readline()
    header4 = lunit.readline()
    header5 = lunit.readline()
    header6 = lunit.readline()
    header7 = lunit.readline()
    header8 = lunit.readline()
    header9 = lunit.readline()
    header10 = lunit.readline()
    header11 = lunit.readline()
    header12 = lunit.readline()
    header13 = lunit.readline()
    header14 = lunit.readline()
    header15 = lunit.readline()

    header16 = lunit.readline()
    header17 = lunit.readline()
    header18 = lunit.readline()
    header19 = lunit.readline()
    header20 = lunit.readline()
    header21 = lunit.readline()
    header22 = lunit.readline()
    header23 = lunit.readline()

    #Lonmin   = int((header10.split())[4])
    #Lonmax   = int((header10.split())[5])
    #Latmin   = int((header11.split())[4])
    #Latmax   = int((header11.split())[5])

    #XSIZE    = int((header12.split())[3])
    #YSIZE    = int((header13.split())[3])

    nbMCS    = int((header17.split())[5])
    data = []
    iMCS = -1
    lines = lunit.readlines()
    
    for iline in lines:
        Values = iline.split()
        if(Values[0] == '==>'):
            #
            # Read the integrated parameters of the convective systems
            ###########################################################
            data.append(MCS_IntParameters())
            iMCS = iMCS+1
            data[iMCS].label 			= int(Values[1])       # Label of the convective system in the segmented images
            data[iMCS].qc_MCS			= int(Values[2])       # Quality control of the convective system 
            data[iMCS].classif			= float(Values[3])		# classification of the MCS (class I/class IIa / Class IIb)
            data[iMCS].nbslot			= float(Values[4])		# duration of the convective system in number of slot

            data[iMCS].duration			= float(Values[4])/2	# duration of the convective system (hr)
            data[iMCS].Utime_Init		= float(Values[5])		# time TU of initiation of the convective system
            data[iMCS].localtime_Init	= float(Values[6])		# local time of inititiation
            data[iMCS].lonInit			= float(Values[7])		# longitude of the center of mass at inititiation
            data[iMCS].latInit			= float(Values[8])		# latitude of the center of mass at inititiation
            data[iMCS].Utime_End		= float(Values[9])		# time TU of dissipation of the convective system
            data[iMCS].localtime_End	= float(Values[10])		# local hour of dissipation
            data[iMCS].lonEnd			= float(Values[11])		# longitude of the center of mass at dissipation
            data[iMCS].latEnd			= float(Values[12])		# latitude of the center of mass at dissipation
            data[iMCS].vavg				= float(Values[13])		# average velocity during its life cycle(m/s)
            data[iMCS].dist				= float(Values[14])		# distance covered by the convective system during its life cycle(km)

            data[iMCS].lonmin			= float(Values[15])		# longitude min of the center of mass during its life cycle
            data[iMCS].lonmax			= float(Values[16])		# longitude max of the center of mass during its life cycle
            data[iMCS].latmin			= float(Values[17])		# latitude min of the center of mass during its life cycle
            data[iMCS].latmax			= float(Values[18])		# latitude max of the center of mass during its life cycle

            data[iMCS].tbmin			= int(Values[19])		# minimum Brigthness temperature (K)

            data[iMCS].surfmaxPix_235K	= int(Values[20])		# maximum surface for a 235K threshold of the convective system during its life cycle (pixel)
            data[iMCS].surfmaxkm2_235K	= float(Values[21])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].surfmaxkm2_220K	= float(Values[22])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].surfmaxkm2_210K	= float(Values[23])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)
            data[iMCS].surfmaxkm2_200K	= float(Values[24])		# maximum surfacefor a 235K threshold of the convective system during its life cycle (km2)

            data[iMCS].surfcumkm2		= float(Values[25]) 	# integrated cumulated surface for a 235K threshold of the convective system during its life cycle (km2)

            data[iMCS].clusters = MCS_Lifecycle()

        else:	
            #
            # Read the parameters of the convective systems 
            # along their life cycles
            ##################################################
            data[iMCS].clusters.qc_im.append(int(Values[0]))       			# quality control on the Infrared image
            data[iMCS].clusters.tbmin.append(int(float(Values[1])))       			# min brightness temperature of the convective system at day TU (K)
            data[iMCS].clusters.tbavg_235.append(int(float(Values[2])))       			#  average brightness temperature of the convective system at day TU (K) 
            data[iMCS].clusters.tbavg_208.append(int(float(Values[3])))
            data[iMCS].clusters.tbavg_200.append(int(float(Values[4])))
            data[iMCS].clusters.tbavg_90th.append(int(float(Values[5])))
            data[iMCS].clusters.Utime.append(float(Values[6]))       			# day TU 
            data[iMCS].clusters.Localtime.append(float(Values[7]))       		# local hour (h)
            data[iMCS].clusters.lon.append(float(Values[8]))       			# longitude of the center of mass (°)
            data[iMCS].clusters.lat.append(float(Values[9]))       			# latitude of the center of mass (°)
            data[iMCS].clusters.x.append(int(float(Values[10])))          			# column of the center of mass (pixel)
            data[iMCS].clusters.y.append(int(float(Values[11])))          			# line of the center of mass(pixel)
            data[iMCS].clusters.velocity.append(float(Values[12]))       		# instantaneous velocity of the center of mass (m/s)

            data[iMCS].clusters.semiminor_220K.append(float(Values[13]))       # semi-minor axis of the ellipse at 220K
            data[iMCS].clusters.semimajor_220K.append(float(Values[14]))       # semi-major axis of the ellipse at 235K
            data[iMCS].clusters.excentricity_220K.append(float(Values[15]))  	# excentricity of the ellipse at 235K
            data[iMCS].clusters.orientation_220K.append(float(Values[16]))   	# orientation of the ellipse at 235K

            data[iMCS].clusters.semiminor_235K.append(float(Values[17]))       # semi-minor axis of the ellipse at 235K
            data[iMCS].clusters.semimajor_235K.append(float(Values[18]))       # semi-major axis of the ellipse at 235K
            data[iMCS].clusters.excentricity_235K.append(float(Values[19]))  	# excentricity of the ellipse at 235K
            data[iMCS].clusters.orientation_235K.append(float(Values[20]))   	# orientation of the ellipse at 235K

            data[iMCS].clusters.surfPix_235K.append(int(float(Values[21])))       		# surface of the convective system at time day TU (pixel)
            data[iMCS].clusters.surfPix_210K.append(int(float(Values[22])))       		# surface of the convective system at time day TU (pixel)

            data[iMCS].clusters.surfkm2_235K.append(float(Values[23]))  		# surface of the convective system for a 235K threshold
            data[iMCS].clusters.surfkm2_220K.append(float(Values[24]))  		# surface of the convective system for a 200K threshold
            data[iMCS].clusters.surfkm2_210K.append(float(Values[25]))  		# surface of the convective system for a 210K threshold
            data[iMCS].clusters.surfkm2_200K.append(float(Values[26]))  		# surface of the convective system for a 220K threshold

    return data   

year = 2016
month = 5
day = 15
    
def parameter_TOOCAN(iregion):
    ##################################################
    ##################################################
    ##################################################
    # DEFINE the Path of the TOOCAN runs 
    ##################################################
    ##################################################
    ##################################################

    if(iregion == 0):
            REGION      = 'AFRICA'
            SAT         = 'MSG3'
            lonmin      = -40
            lonmax      = 40
            latmin      = -20
            latmax      = 20
    if(iregion == 1):
            REGION      = 'INDIA'
            SAT         = 'MET7'
            lonmin      = 40
            lonmax      = 101
            latmin      = -20
            latmax      = 20
    #
    # if Region == WESTERNPACIFIC and the date is before 20150601
    # then the satellite used is MTSAT-2
    # For technical reasons inherent to the Geostationary satellite, the SouthScan is not processed
    ###############################################################################################
    if( (iregion == 2) & (year*100+month < 201506)) :
            REGION      = 'WESTERNPACIFIC'
            SAT         = 'MTSAT2'
            lonmin      = 101
            lonmax      = 185
            latmin      = 0
            latmax      = 20
    #
    # if Region == WESTERNPACIFIC and the date is after 20150601
    # then the satellite used is HIMAWARI-8
    ###############################################################
    if( (iregion == 2) & (year*100+month >= 201506)) :
            REGION      = 'WESTERNPACIFIC'
            SAT         = 'HIMAWARI8'
            lonmin      = 101
            lonmax      = 185
            latmin      = 0
            latmax      = 20		
    if(iregion == 3):
            REGION      = 'EASTERNPACIFIC'
            SAT         = 'GOES15'
            lonmin      = 185
            lonmax      = 250
            latmin      = -20
            latmax      = 20
    if(iregion == 4):
            REGION      = 'AMERICA'
            SAT         = 'GOES13'
            lonmin      = -106
            lonmax      = -40
            latmin      = -20
            latmax      = 20
            
    return REGION, latmin, latmax, lonmin, lonmax      



