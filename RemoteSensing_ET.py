# -*- coding: utf-8 -*-
#author: curtis crandall 2017

# Imports
import pandas as pd
import numpy as np
import arcpy
import os


# File import; change the names depending on which image to process
currentDate = 20170527
dayInYear = 147
img_csv = '20170527_EVI_CALD.csv'  # imported csv with evi values as 'grid_code'
img_csv2 = "20170527_ET_CALD.csv"  # new csv with calculated et
img_csv3 = '20170527_NDVI_CALD.csv'  # imported csv with ndvi values
img_csv4 = "Agrimet_CALD_modelData_current.csv"
img_csv5 = '20170527_Albedo_CALD.csv'  # csv with albedo values
img_csv10 = "Agrimet_CALD_modelData_appends.csv"
csvRowCount = sum(1 for line in open('20170527_EVI_CALD.csv'))  # Get row count for array sizing


# Open a data frame of evi values using pandas
df = pd.read_csv(img_csv, names=['OBJECTID', 'pointid', 'EVI'], skiprows=1)
EVIs = df['EVI'].values
# Open a data frame of ndvi values using pandas
df2 = pd.read_csv(img_csv3, names=['OBJECTID', 'pointid', 'NDVI'], skiprows=1)
dfNDVIs = df2['NDVI'].values
# Open a data frame of albedo
df3 = pd.read_csv(img_csv5, names=['OBJECTID', 'pointid', 'EVI', 'ALBEDO'], skiprows=1)
ALBEDOs_grid = df3['ALBEDO'].values

# Agrimet Station Inputs
# AgMet_ETkm_input = 0  # k.m. reference ET inches per day
AgMet_ETrs_input = 0.26  # ET alfalfa reference value inches per day
AgMet_ETos_input = 0.22  # ET grass reference value inches per day
Tmin_input = 45.64  # (MN) input temperature data in °F
Tmax_input = 81.70  # (MX) input temperature data in °F
Tmean_input = 63.90  # (MM) input temperature data in °F
Tmean3_input = 60.75  # (MM) previous 3 day average temp in F
Rnet_input = 743.32  # (SR) input solar radiation
RH_mean_input = 51.2  # (TA) mean daily relative humidity as %
Tdew_input = 43.02  # (YM) input dew point temperature data in °F
uz_input = 2.7  # (UA) measured wind speed 2m above the ground surface, mph;
z = 825  # elevation of station in m

# Input Conversions
# AgMet_ETkm = AgMet_ETkm_input*25.4  # kimberly penman-monteith et in mm/day from Agrimet station
AgMet_ETrs = AgMet_ETrs_input*25.4  # ASCE alfalfa reference et in mm/day from Agrimet station
AgMet_ETos = AgMet_ETos_input*25.4  # ASCE grass reference et in mm/day from Agrimet station
uz = uz_input*0.44704  # daily average wind speed m/s
Tmax = (Tmax_input-32)/1.8  # max daily air temperature, °C
Tmin = (Tmin_input-32)/1.8  # min daily air temperature, °C
Tmean = (Tmean_input-32)/1.8  # mean daily air temperature, °C
Tmean3 = (Tmean3_input-32)/1.8  # mean daily air temperature, °C
Tdew = (Tdew_input-32)/1.8  # dew point temperature, °C

# Other Variables
es = 0.611*np.exp((17.27*Tmean)/(Tmean + 237.3))  # saturation vapor pressure, kPa
ea = 0.611*np.exp((17.27*Tdew)/(Tdew + 237.3))  # preferred actual vapor pressure calculation, kPa
svp_def = es-ea  # saturation vapor pressure deficit, kPa
bc = (4.9*10**(-9))*((Tmean)**4)  # stephan boltzman constant
P = 101.3*((293 - (.0065*z))/293)**5.26  # good; atmospheric pressure at elevation z
gamma = .000665*P  # good; psychrometric constant, kPa °C-1
h = 2  # height of the measurement above the ground surface, m.
u2 = uz*(4.87/(np.log(67.8*h - 5.42)))  # wind speed at 2 m height, m s-1
delta = (2503*np.exp((17.27*Tmean)/(Tmean + 237.3)))/((Tmean + 237.3)**2)  # slope of vapor pressure curve, kPa/ºC-1
DT = delta/(delta+gamma*(1+.34*u2))  # Delta Term used to calculate Radiation Term
PT = gamma/(delta+gamma*(1+.34*u2))  # Psi Term used to calculate Wind Term
TT = (900/(Tmean+273))*u2  # Temperature Term used to calculate Wind Term


# NDVI values; for now limiting values of zero to avoid nan ET values
NDVIs = []
for i in range(len(EVIs)):
    x = dfNDVIs[i]
    if x <= 0:
        x = 0
    NDVIs.append(x)
NDVIs = map(float, NDVIs)


# Albedo values; for now limiting values of zero to avoid nan ET values
ALBEDOs = []
for i in range(len(EVIs)):
    x = ALBEDOs_grid[i]
    if x <= .01:
        x = .01
    ALBEDOs.append(x)
ALBEDOs = map(float, ALBEDOs)


# Leaf Area Index values; for now limiting values of zero to avoid nan ET values
LAIs = []
for i in range(len(EVIs)):
    x = 10.*(NDVIs[i]**3.5)
    if x > 6.:
        x = 6.
    if x <= 0:
        x = 0.1
    LAIs.append(x)
LAIs = map(float, LAIs)


# Net Solar Radiation values; for now limiting values of zero to avoid nan ET values
RNs = []
for i in range(len(EVIs)):
    x = Rnet_input*0.04184*(1.-ALBEDOs[i])
    if x <= 0:
        x = 0.1
    RNs.append(x)
RNs = map(float, RNs)


# Calculate Soil Heat Flux; for now limiting values of zero to avoid nan ET values
HEATFLUXs = []
for i in range(len(EVIs)):
    x = 0.1 + 0.17 * math.exp((-.55)*LAIs[i])
    if x < 0:
        x = 0
    HEATFLUXs.append(x)
HEATFLUXs = map(float, HEATFLUXs)


# Alfalfa Reference ET
ETrs = []
for i in range(len(EVIs)):
    x = (0.408*delta*(RNs[i]-HEATFLUXs[i]) + gamma*(1600./(Tmean+273.))*u2*svp_def)/(delta + gamma*(1.+.38*u2))
    if x < 0:
        x = 0
    ETrs.append(x)
#ETrs = map(float, ETrs)

# Calculate Actual ET from Alfalfa Reference ET
EVI_ETra = []
for i in range(len(EVIs)):
    c = 2.5  # fitting coefficient to 20% of reported alfalfa reference for this pixel
    x = ETrs[i]*1.65*(1-math.exp((-2.25)*EVIs[i]))-c
    if x < 0:
        x = 0
    EVI_ETra.append(x)
#EVI_ETra = map(float, EVI_ETra)

# Grass Reference ET
ETos = []
for i in range(len(EVIs)):
    x = (0.408*delta*(RNs[i]-HEATFLUXs[i]) + gamma*(900./(Tmean+273.))*u2*svp_def)/(delta + gamma*(1.+.34*u2))
    if x < 0:
        x = 0
    ETos.append(x)
#ETos = map(float, ETos)

# Calculate Actual ET from Grass Reference ET
EVI_EToa = []
for i in range(len(EVIs)):
    c = 2.5  # fitting coefficient
    x = ETos[i]*1.65*(1-math.exp((-2.25)*EVIs[i]))-c
    if x < 0:
        x = 0
    EVI_EToa.append(x)
#EVI_EToa = map(float, EVI_EToa)

# Calculate ETrF; the current ET rate expressed as a fraction of maximum reference ET
# ETrF = calculated actual / reference
# Alfalfa ETrF
ETrF_rs = []
for i in range(len(EVIs)):
    x = abs((EVI_ETra[i])/(ETrs[i]))
    ETrF_rs.append(x)
ETrF_rs = map(float, ETrF_rs)
# Grass ETrF
ETrF_os = []
for i in range(len(EVIs)):
    x = abs((EVI_EToa[i])/(ETos[i]))
    ETrF_os.append(x)
ETrF_os = map(float, ETrF_os)

# Create and populate a new csv to import back into arcmap
df4 = pd.DataFrame({"ETrF_rs": ETrF_rs, "ETrF_os": ETrF_os, "ETra": EVI_ETra, "EToa": EVI_EToa, "EVI": EVIs,
                    "NDVI": NDVIs, "RN": RNs, "LAI": LAIs, "Albedo": ALBEDOs, "SoilHeatFlux": HEATFLUXs})
df4.to_csv(img_csv2, index=False)

# BOII Agrimet specific pixel data; nutty work-around to extract attributes of an index
BOII_Agrimet = df4.ix[58772]
BOII_albedo = []
BOII_EToa = []
BOII_ETrF_os = []
BOII_ETrF_rs = []
BOII_ETra = []
BOII_EVI = []
BOII_LAI = []
BOII_NDVI = []
BOII_RN = []
BOII_SoilHeatFlux = []
for i in range(len(BOII_Agrimet)-9):
    a = BOII_Agrimet[0]
    BOII_albedo.append(a)
    b = BOII_Agrimet[1]
    BOII_EToa.append(b)
    c = BOII_Agrimet[2]
    BOII_ETrF_os.append(c)
    d = BOII_Agrimet[3]
    BOII_ETrF_rs.append(d)
    e = BOII_Agrimet[4]
    BOII_ETra.append(e)
    f = BOII_Agrimet[5]
    BOII_EVI.append(f)
    g = BOII_Agrimet[6]
    BOII_LAI.append(g)
    h = BOII_Agrimet[7]
    BOII_NDVI.append(h)
    j = BOII_Agrimet[8]
    BOII_RN.append(j)
    k = BOII_Agrimet[9]
    BOII_SoilHeatFlux.append(k)

df5 = pd.DataFrame({"Date": currentDate, "_ETrF_rs": BOII_ETrF_rs, "_ETrF_os": BOII_ETrF_os, "_ETra": BOII_ETra,
                    "_EToa": BOII_EToa, "_EVI": BOII_EVI, "_NDVI": BOII_NDVI, "_RN": BOII_RN, "_LAI": BOII_LAI,
                    "_Albedo": BOII_albedo, "_SoilHeatFlux": BOII_SoilHeatFlux})
df5.to_csv(img_csv4, index=False)
# Append current station output to list
# Only run this one time per image date
with open(img_csv10, 'a') as f:
    df5.to_csv(f, header=False, index=False)


# print df5 to check pixel stats
print df4.ix[65564]
print df4.ix[58772]

