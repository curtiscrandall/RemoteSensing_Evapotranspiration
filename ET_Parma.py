# -*- coding: utf-8 -*-
"""
Created on Mon Aug 06 12:29:07 2018

@author: crc
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 21 10:31:33 2018

@author: crc
"""

# imports for script
import numpy as np


def etProcessing(station, year, month, day):
    
    # imports for function
    from osgeo import gdal
    from osgeo.gdalconst import GDT_Float32
    import numpy as np
    import sys
    import pandas as pd


    # -----------------------------------------------------------------------------------------------------
    # IMAGE FILE LOCATIONS
    # -----------------------------------------------------------------------------------------------------
    # ndviRaster   = "D:\GIS_Data\Thesis\Scripts\Rasters\ndvi0623.tif"
    # albedoRaster = "D:\GIS_Data\Thesis\Scripts\Rasters\albedo0623.tif"
    eviRaster  = "D:\GIS_Data\Thesis\Scripts\Rasters\evi{}{}.tIF" .format(month, day)
    newFile    = "D:\GIS_Data\Thesis\Scripts\Rasters\etAg{}{}.tif" .format(month, day)
    
    # Register files for all of the GDAL drivers
    gdal.AllRegister()
    
    # Print Date
    # print('----Starting: {}{}{}{}----' .format(station, year, month, day))
    
    # -----------------------------------------------------------------------------------------------------
    # CLIMATE DATA FROM AGRIMET STATIONS
    # -----------------------------------------------------------------------------------------------------
    
    # Parma Agrimet station information
    if station == "Parma":
        # print("--PARMA INFORMATION--")
        url = "https://www.usbr.gov/pn-bin/daily.pl?station=PMAI&year=%s&month=%s&day=%s&year=%s&month=%s&day=%s&pcode=ETRS&pcode=MN&pcode=MX&pcode=MM&pcode=SR&pcode=TA&pcode=YM&pcode=UA" % (year, month, day, year, month, day)
        # table from website
        urlCSV         = pd.read_table(url)
        urlValues      = urlCSV[18:20].values
        AgMet_KM_input = float(urlValues[1,0][19:23])
        Tmin_input     = float(urlValues[1,0][31:36])
        Tmax_input     = float(urlValues[1,0][43:49])
        Tmean_input    = float(urlValues[1,0][57:62])
        Rnet_input     = float(urlValues[1,0][69:75])
        RH_mean_input  = float(urlValues[1,0][83:88])
        Tdew_input     = float(urlValues[1,0][96:101])
        uz_input       = float(urlValues[1,0][110:114])
        z              = 703

    # -----------------------------------------------------------------------------------------------------
    # EVI RASTER ARRAY
    # -----------------------------------------------------------------------------------------------------
    # open the image
    inDs = gdal.Open(eviRaster)
    if inDs is None:
      print 'Could not open image file'
      sys.exit(1)
    
    # read in data for first band
    band1 = inDs.GetRasterBand(1)
    # calculate rows
    rows = inDs.RasterYSize
    # calculate columns
    cols = inDs.RasterXSize
    # create array
    eviData = band1.ReadAsArray(0, 0, cols, rows)
    
    # Use identified pixel points from Arcmap of sites
    # Array coordinates of weather stations within EVI raster
    eviParma = eviData[5034,537]
    # print('PMAI station: %s' %eviParma)
    # eviBoise = eviData[7220,6645]
    # print('BOII station: %s' %eviBoise)
    
    
    # -----------------------------------------------------------------------------------------------------
    # ET CALCULATIONS
    # -----------------------------------------------------------------------------------------------------
    # Input Conversions
    
    AgMet_KM = AgMet_KM_input*25.4  # ASCE alfalfa reference et in mm/day from Agrimet station
    uz = uz_input*0.44704  # daily average wind speed m/s
    Tmax = (Tmax_input-32)*.5556  # max daily air temperature, °C
    Tmin = (Tmin_input-32)*.5556  # min daily air temperature, °C
    Tmean = (Tmean_input-32)*.5556  # mean daily air temperature, °C
    Tmean = (Tmin + Tmax)/2  # mean daily air temperature, °C
    Tdew = (Tdew_input-32)*.5556  # dew point temperature, °C
    
    # Other Variables
    k       = 273.15
    G       = .1  # soil heat flux; check metric for local derivation
    a       = .25  # albedo
    # es      = 0.6108*np.exp((17.27*Tmean)/(Tmean + 237.3))  # saturation vapor pressure, kPa
    es      = (0.6108*np.exp((17.27*Tmin)/(Tmin + 237.3)) + 0.6108*np.exp((17.27*Tmax)/(Tmax + 237.3)))/2  # saturation vapor pressure, kPa
    ea      = 0.6108*np.exp((17.27*Tdew)/(Tdew + 237.3))  # preferred actual vapor pressure calculation, kPa
    svp_def = es-ea  # saturation vapor pressure deficit, kPa
    bc      = (5.67*(10**(-8)))/((k)**4)  # stephan boltzman constant
    P       = 101.3*((293 - (.0065*z))/293)**5.26  # atmospheric pressure at elevation z
    gamma   = .000665*P  # psychrometric constant, kPa °C-1
    h       = 2  # height of the measurement above the ground surface, m.
    u2      = uz*(4.87/(np.log(67.8*h - 5.42)))  # wind speed at 2 m height, m s-1
    delta   = (2503*np.exp((17.27*Tmean)/(Tmean + 237.3)))/((Tmean + 237.3)**2)  # slope of vapor pressure curve, kPa/ºC-1
    # DT      = delta/(delta+gamma*(1+.34*u2))  # Delta Term used to calculate Radiation Term
    # PT      = gamma/(delta+gamma*(1+.34*u2))  # Psi Term used to calculate Wind Term
    # TT      = (900/(Tmean+273))*u2  # Temperature Term used to calculate Wind Term
    # Rns      = Rnet_input*0.04184*(1. - a)
    # Rnl     = bc*0.2*(0.34-0.14*np.sqrt(ea))*(((Tmax + k)**4+(Tmin + k)**4)/2)
    # Lin     = (2.7*ea + .245*(Tmean+k)-45.14)
    # Lout    = 0.98*bc*((Tmean+k)**4)+(1-0.98)*Lin
    # Rnl     = 0.98*bc*(Tmean+k)**4
    # Rnet    = Rns + Rnl
    rnS2    = .7  # coefficient to convert global rad to net rad (Alados et al. 2003)
    rn      = rnS2 * Rnet_input * 0.04184 * (1. - a)
    Rnet    = rnS2 * (Rnet_input - 48.5) * 0.04184 * (1. - a)
    
    
    # Crop Coefficients for cool season grass/turf (kentucky blue grass)
    if eviParma < .5:
        cc = .3
    if eviParma >= .5:
        cc = .95
    #if eviParma >= .7:
    #    cc = .8
    #if eviParma >= .8:
    #    cc = 1
    # Dual Crop Coefficient
    grassHeight = 1.5  # meters
    rh = (ea/es)  # relative humidity
    dc = cc + (0.04*(uz-2)-0.004*(rh-45))*((grassHeight/3)**0.3)
    
    # ET Alfalfa Standardized Reference
    ETrs = (0.408*delta*(Rnet - G) + gamma*(1600./(Tmean+273.))*u2*svp_def)/(delta + gamma*(1.+.38*u2))
    # Difference between calculated and reported reference ET
    dif = abs(1- ETrs / AgMet_KM)*100
    dif2 = abs(ETrs - AgMet_KM)
    # ET Actual using crop coefficient
    # ETcc = cc*AgMet_KM
    ETdc = dc*AgMet_KM
    # print('Agmet Ref: %.2f' % AgMet_KM)
    # print('Model Ref: %.2f' % ETrs)
    # print('Ref ET Difference: %.2f percent' % dif)
    
    # ET ACTUAL FROM EVI
    # Calibration
    pixel = eviParma
    c = 0.49  # fitting coefficient
    b = 1.48  # fitting coefficient 1.48
    ETparma = ETrs * b * (1. - np.exp((-2.25) * pixel)) - c
    # Correct to match dual crop coefficient
    # c = ETparma - ETdc  # fitting coefficient
    # ETparma = ETrs * b * (1. - np.exp((-2.25) * pixel)) - c
    
    # Difference of Eta and Model ETa
    # print('Model Alfalfa Reference ET: %.2f' % ETrs)
    # print('Grass ETa: %.2f' % ETdc)
    print('%.2f,%.2f' % (ETdc,ETparma))  # just print ref and act et for station for excel import
    # print('Old Rnet: %.2f' %rn)
    # print('New Rnet: %.2f' %Rnet)
    # print('cc Model ETa: %.2f' % ETboiicc)
    # print('BOII Avg Model ETa: %.2f' % ETboiiavg)
    # print('PMAI Avg Model ETa: %.2f' % ETpmaiavg)
    # print('Model calibrated coefficient: %.2f' % c)
    # diff3 = abs(1- ETboiicc / ETboiiavg)*100
    # print('Difference of ETa and Model ET: %.2f' % diff3)

    """
    # -----------------------------------------------------------------------------------------------------
    # CREATE NEW ET IMAGE
    # -----------------------------------------------------------------------------------------------------
    # create the output image
    driver = inDs.GetDriver()
    # print driver
    outDs = driver.Create(newFile, cols, rows, 1, GDT_Float32)
    # Process whole raster
    outNewData = ETrs * b * (1. - np.exp((-2.95) * eviData)) - c
    outNewData[outNewData <= 0.01] = 0.01  # originiall = -1
    
    # write the data; assign nodata values
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(outNewData)
    outBand.SetNoDataValue(-1)
    # georeference the image and set the projection
    outDs.SetGeoTransform(inDs.GetGeoTransform())
    outDs.SetProjection(inDs.GetProjection())
    outDs = None    
    # delete current stored mem data
    del eviData, outDs, outBand, outNewData
    """
    # print('Finished: {}{}{}' .format(year, month, day))
    return dif, dif2, c, ETparma, AgMet_KM, station


# Processing
site       = "Parma"
etOutput   = []
# Process Each FULL Image using the et function
et03 = etProcessing(site, "2017", "05", "27")
etOutput.append(et03)
et05 = etProcessing(site, "2017", "06", "26")
etOutput.append(et05)
et07 = etProcessing(site, "2017", "07", "06")
etOutput.append(et07)
et09 = etProcessing(site, "2017", "07", "16")
etOutput.append(et09)
et11 = etProcessing(site, "2017", "07", "26")
etOutput.append(et11)
et13 = etProcessing(site, "2017", "08", "05")
etOutput.append(et13)
et15 = etProcessing(site, "2017", "08", "15")
etOutput.append(et15)
et17 = etProcessing(site, "2017", "08", "25")
etOutput.append(et17)
et20 = etProcessing(site, "2017", "09", "24")
etOutput.append(et20)
et21 = etProcessing(site, "2017", "10", "14")
etOutput.append(et21)
et22 = etProcessing(site, "2017", "10", "24")
etOutput.append(et22)


import pandas as pd
df = pd.DataFrame(etOutput, columns=["ETr/Model", "Model-ETr", "coef", "ETpmai", "AgMet_KM", "station"])
df.to_csv('D:\GIS_Data\Thesis\Scripts\Rasters\etOutput{}.csv' .format(site), index=False)


# Simple stats using the et function outputs of each image
etDiff    = [et03[0], et05[0], et07[0], et09[0], et11[0], et13[0], et15[0], et17[0], et20[0], et21[0], et22[0]]
etSTD     = [et03[1], et05[1], et07[1], et09[1], et11[1], et13[1], et15[1], et17[1], et20[1], et21[1], et22[1]]
etC       = [et03[2], et05[2], et07[2], et09[2], et11[2], et13[2], et15[2], et17[2], et20[2], et21[2], et22[2]]
etEviSTD  = [et03[4], et05[4], et07[4], et09[4], et11[4], et13[4], et15[4], et17[4], et20[4], et21[4], et22[4]]
etDiff    = np.array(etDiff)
etSTD     = np.array(etSTD)
etC       = np.array(etC)
print('----Mean Modeling Coefficient----')
print(np.mean(etC))
# print('----Mean Difference in Ref ET----')
# print(np.mean(etDiff))
print('----Standard Deviation of ET----')
print(np.std(etEviSTD))






