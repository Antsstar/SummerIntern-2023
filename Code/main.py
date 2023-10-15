import numpy as np
import math as m
import pandas as pd
import SatelliteClass as sc
import TerminalClass as tc
import ChannelClass as cc
import VectorClass as vc #new

import scipy.constants as spc
from pathlib import Path

from skyfield.api import load, wgs84 #new

def LinkBudget(destination, origin, bandwidth, elevation):
    k_boltz = -228.6
    testCC = cc.Channel(destination, origin, bandwidth, elevation)
    SNR = testCC.EIRP - testCC.BandWidth - testCC.get_FSPL() - testCC.get_AtmosLoss() + testCC.G_T - k_boltz
    return SNR


####Testing Integration
stations_url = 'OneWebDatabase.txt'
satellites = load.tle_file(stations_url)
print('Loaded', len(satellites), 'satellites')

ts = load.timescale()

cities={'London':[51.5072, -0.1276],'Birmingham':[52.4862, -1.8904]}
#cities={'London':[51.5072, -0.1276],'Birmingham':[52.4862, -1.8904], 'Belfast': [54.5973, -5.9301], 'Cardiff': [51.4837, -3.1681], 'Edinburgh': [55.9533, -3.1883], 'Warrington': [53.3900, -2.5970], 'Baldock': [51.9895, -0.1891], 'Manchester': [53.4808, -2.2426]}

Locations = {}
for city,latlon in zip(cities.keys(),cities.values()):
    Locations[city]=wgs84.latlon(*latlon)

#print(Locations)

something = {} #Dictionary of all Satellites with a Skyfield and SatClass pair of instances.
for sat in satellites:
    something[sat.name] = {"Skyfield": sat, "SatClass": sc.Satellite(sat.name, downlink=11.7, uplink=14.25, eirp=-13.4, g_t=-1)}

## Simulation period. Range function contains minutes
t = ts.utc(2023, 7, 24, 0, range(0, 61)) #24*60

column_list = ["Time_Index", "Time_in_UTC", "Ground_Station", "GS_Vector", "Serving_Satellite", "SS_Vector", "GS-SS_vector", "Elevation_Angle", "Distance_from_GS_(km)", "SNR_DL", "SNR_UL", "Visible_Sats"]
df = pd.DataFrame(columns=column_list)

for index_ti,ti in enumerate(t):
    if index_ti not in [0, 10, 20, 30, 40, 50, 60]:
        continue

    for sat_name, classes in zip(something.keys(), something.values()):
        lat, long = wgs84.latlon_of(classes["Skyfield"].at(ti))
        elev = wgs84.height_of(classes["Skyfield"].at(ti))
        geo = [long.degrees, lat.degrees, elev.km]
        something[sat_name]["SatClass"].set_Position(classes["Skyfield"].at(ti).position.km, geo)

    Visible = {}
    for city, place in zip(Locations.keys(), Locations.values()):
        Visible[city] = {sat.name: [something[sat.name]["SatClass"].get_Position(), (sat-place).at(ti).altaz()[-1].km,(sat-place).at(ti).altaz()[0].degrees] for sat in satellites if (sat - place).at(ti).altaz()[0].degrees > 25}

    max_SNRs = {}

    for city, options in zip(Visible.keys(), Visible.values()):
        testTC = tc.Terminal(city, cities[city])
        testTC.set_Position(Locations[city].at(ti).position.km)
        SNRs = {}
        for sat_name in options.keys():
            #print(options[sat_name][-1])
            SNRs[sat_name] = LinkBudget(testTC, something[sat_name]["SatClass"], 4e3, options[sat_name][-1])
            #break
        #print(SNRs)
        max_SNRs[city] = {max(SNRs, key=SNRs.get): max(SNRs.values())}

    for city in cities.keys():
        GS_Vect = vc.Vector(*(Locations[city].at(ti).position.km), 'cartesian')
        Ser_Sat = [dummy for dummy in max_SNRs[city].keys()]
        Ser_Sat = Ser_Sat[0]
        SS_Vect = something[Ser_Sat]["SatClass"].get_Position()
        El_angle = Visible[city][Ser_Sat][-1]
        Dist_SS = Visible[city][Ser_Sat][1]
        SS_Snr = [dummy for dummy in max_SNRs[city].values()]
        SS_Snr = SS_Snr[0]
        testTC = tc.Terminal(city, cities[city])
        testTC.set_Position(Locations[city].at(ti).position.km)
        Snr_UL = LinkBudget(something[Ser_Sat]["SatClass"], testTC, 125e6, El_angle)

        df.loc[df.shape[0]+1] = [index_ti, ti, city, GS_Vect, Ser_Sat, SS_Vect, GS_Vect - SS_Vect, El_angle, Dist_SS, SS_Snr, Snr_UL, {sat_name: {"Earth": something[sat_name]["SatClass"].get_Position() ,"GroundSat": GS_Vect - something[sat_name]["SatClass"].get_Position()} for sat_name in Visible[city].keys()} ]
        
    
DL_bandwidth=250e6
UL_bandwidth=125e6
mcs=pd.read_csv(Path('DVBS2-table-meta.csv'))

sE_list_DL=[]
modulation_list_DL=[]
sE_list_UL=[]
modulation_list_UL=[]
for i in df.index:
    # DL
    sE_list_DL.append(mcs[mcs.SNR<df.loc[i,'SNR_DL']].SE.max())
    modulation_list_DL.append(mcs.loc[mcs[mcs.SNR<df.loc[i,'SNR_DL']].SE.argmax(),'code'])
    #UL
    sE_list_UL.append(mcs[mcs.SNR<df.loc[i,'SNR_UL']].SE.max())
    modulation_list_UL.append(mcs.loc[mcs[mcs.SNR<df.loc[i,'SNR_UL']].SE.argmax(),'code'])

df['DL_capacity']=np.array(sE_list_DL)*DL_bandwidth/1e6
df['DL_Mode']=modulation_list_DL
df['UL_capacity']=np.array(sE_list_UL)*UL_bandwidth/1e6
df['UL_Mode']=modulation_list_UL

df.to_pickle("Results.pkl") #Saves

