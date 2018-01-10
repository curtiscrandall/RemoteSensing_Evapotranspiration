# -*- coding: utf-8 -*-
#author: curtis crandall
"""
This code takes csv file input, from satellite imagery processed in Arcmap to produce EVI values, to compute EVI based
evapotranspiration over the point associated data from Arcmap. It will produce a new csv that can be imported back into
arcmap and Joined with the original point file based on OBJECTID. This will associated each point (30m pixel) with a
calculated average evapotranspiration value for that day. The numbers can then be correlated with Agrimet station data
in the same pixel.
"""
import csv, math, pandas as pd, numpy as np, matplotlib.mlab as mlab, matplotlib.pyplot as plt

#  use this to extract evi variables from processed image point file csv (arcmap raster to point output)
img_csv = 'L8_20130801_EVI_BOII_2km.csv'  # file name of imported csv with evi values
img_csv2 = "L8_20130801_EVI_BOII_2km_ET.csv"  # file name of created csv with evi values and calculated et
EVIs = []  # list to be populated with grid_code values from the arcmap point file csv

with open(img_csv) as csv_in:  # populates EVIs list with evi values from arcmap csv
    readcsv = csv.reader(csv_in, delimiter=',')  # reads in import csv with evi values
    next(readcsv, None)  # skip header row in csv
    for column in readcsv:  # for columns in the csv
        evi = column[2]  # takes every value from column 3 in the csv into variable evi
        EVIs.append(evi)  # append evi variables from column 3 to EVIs list

# Evapotranspiration computations
# inputs from weather station
"""
SAMPLE DATA
DATE      , ET K.P(in), ASCE Alfalfa(in), ASCE Grass(in), MinTemp(F), MaxTemp(F), MeanTemp(F), SolarRad(langleys), MeanHumid(%), MeanDew(F), MeanWind(mph)
08/01/2013, 0.26,       0.23,             0.20,           55.52,      90.10,      74.80,       623.00,             40.41,        46.27,      1.54
"""
evi = .35  # sample EVI calculated from Landsat 8 on Aug 1, 2013 from the Agrimet station
ET_kp_ref_input = 0.26
ET_alf_ref_input = 0.23
ET_gra_ref_input = 0.20
Rnet_input = 623  # input solar radiation
RH_mean_input = 40.41  # mean daily relative humidity as %
Tmax_input = 90.1  # input temperature data in °F
Tmin_input = 55.52  # input temperature data in °F
Tmean_input = 74.80  # input temperature data in °F
Tdew_input = 46.27  # input dew point temperature data in °F
uz_input = 1.54  # measured wind speed 2m above the ground surface, mph;
z = 828  # elevation of station in m
G = Rnet_input*0.04184*0.1  # soil heat flux density, MJ m-2 d-1
albedo = 0.23  # standardized constant for average albedo
ET_kp_ref = ET_kp_ref_input*25.4  # kimberly penman et in mm/day from agrimet station
ET_alf_ref = ET_alf_ref_input*25.4  # ASCE alfalfa reference et in mm/day from agrimet station
ET_gra_ref = ET_gra_ref_input*25.4  # ASCE grass reference et in mm/day from agrimet station
uz = uz_input*0.44704  # daily average wind speed m/s
Tmax = (Tmax_input-32)/1.8  # max daily air temperature, °C
Tmin = (Tmin_input-32)/1.8  # min daily air temperature, °C
Tmean = (Tmean_input-32)/1.8  # mean daily air temperature, °C
Tdew = (Tdew_input-32)/1.8  # dew point temperature, °C
es = 0.611*math.exp((17.27*Tmean)/(Tmean + 237.3))  # saturation vapor pressure, kPa
ea = 0.611*math.exp((17.27*Tdew)/(Tdew + 237.3))  # preferred actual vapor pressure calculation, kPa
svp_def = es-ea  # saturation vapor pressure deficit, kPa
bc = (4.9*10**(-9))*((Tmean)**4)  # stephan boltzman constant
# R_Lin = .98*bc*(Tmean+273.15)**4
# R_Lin2 = 2.7*ea+.245*(Tmean + 273.15)-45.14
Rnet = Rnet_input*0.04184*(1-albedo)  # net radiation at the crop surface, langleys to MJ m-2 d-1
P = 101.3*((293 - (.0065*z))/293)**5.26  # good; atmospheric pressure at elevation z
gamma = .000665*P  # good; psychrometric constant, kPa °C-1
h = 2  # height of the measurement above the ground surface, m.
u2 = uz*(4.87/(math.log(67.8*h - 5.42)))  # wind speed at 2 m height, m s-1
delta = (2503*math.exp((17.27*Tmean)/(Tmean + 237.3)))/((Tmean + 237.3)**2)  # slope of vapor pressure curve, kPa/ºC-1
DT = delta/(delta+gamma*(1+.34*u2))  # Delta Term used to calculate Radiation Term
PT = gamma/(delta+gamma*(1+.34*u2))  # Psi Term used to calculate Wind Term
TT = (900/(Tmean+273))*u2  # Temperature Term used to calculate Wind Term
ETos = (0.408*delta*(Rnet-G) + gamma*(900/(Tmean+273))*u2*svp_def)/(delta + gamma*(1+.34*u2))  # reference et, mm/day
ETevi = ETos*1.65*(1-math.exp(-2.25*evi))-0.169  # calculate ET based on EVI values in same pixel
Difference1 = (1-ET_gra_ref/ETos)*100  # percentage difference between my FAO ET and Agrimet grass ref ET
Difference2 = (1-ET_gra_ref/ETevi)*100 # percentage difference between my EVI ET and Agrimet grass ref ET

ET_evi = []  # list with calculated ET values using EVI
for i in range(len(EVIs)):  # calculate et values for each evi value in the EVIs list
    x = ETos*1.65*(1-math.exp(-2.25*float(EVIs[i])))-0.169  # ET calculation using EVI values
    ET_evi.append(x)  # append evi calculated et to ET_evi list
df = pd.DataFrame({"ET_mm": ET_evi, "EVI": EVIs})  # Create csv from EVIs and Calculated ET
df.to_csv(img_csv2, index=False)  # create new csv with EVI and calculated EVI ET

# create histogram from ET points
# numpy.histogram(a, bins=10, range=None, normed=False, weights=None, density=None)
# ET_hist = np.histogram(ET_evi, bins=30, range=None, normed=False, weights=None, density=None)

# print list
print ('Agrimet ET for grass is %.2f mm/day' % float(ET_gra_ref))
print ('Calculated ET using EVI value is %.2f ' % float(ETevi))
print ('Difference is %.1f percent' % float(Difference2))



