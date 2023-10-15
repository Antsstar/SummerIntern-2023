import numpy as np
import math as m
import pandas as pd
import VectorClass as vc
import SatelliteClass as sc
import TerminalClass as tc

import itur
import astropy.units as u


class Channel: #Channel class for storing properties of the active connection channel.
    """ Define an Active Connection
    Key Attributes:
    Beam direction - Direction and magnitude of signal propagation defined in VectorClass (km)
    Elevation Angle - Stores the elevation angle respect to the Earth's Horizon of the Beam trajectory (degrees)
    BandWidth - Stores the Bandwidth frequency of the active connection in (dBHz)

    Key Functions:
    __Calculate_FSPL - Finds the Free-Space Path Loss (dB)
    __Calculate_AtmosLoss - Finds the Atmospheric Loss using ITU-R model (dB)
 """
    __C = 3 * 10**(-4)

    def __init__(self, Destination, Origin, BandWidth, Elevation):
        self.__Beam_Direction = (Destination.get_Position() - Origin.get_Position())
        self.__Start_Position = Origin.get_Position()
        self.__Elevation_Angle = Elevation
        self.BandWidth = 10*m.log10(BandWidth)
        self.__angle = Origin.get_Pointing().angle(self.__Beam_Direction)
        self.EIRP = Origin.Properties(self.__angle)[-1]
        self.G_T, self.G = Destination.Properties(self.__angle)[:-1]

        if isinstance(Origin, sc.Satellite):
            self.__Base_Frequency = Origin.get_Downlink()
            self.Sat_Geo = Origin.get_Geo()
            self.GS_Geo = Destination.get_Geo()
        else:
            self.__Base_Frequency = Destination.get_Uplink()
            self.Sat_Geo = Destination.get_Geo()
            self.GS_Geo = Origin.get_Geo()
        
        self.__Calculate_FSPL()
        self.__Calculate_AtmosLoss()

    def __Calculate_FSPL(self):
        self.__FSPL = 20*m.log10(self.__Beam_Direction.magnitude()) + 20*m.log10(self.__Base_Frequency) + 20*m.log10((4*m.pi) / self.__C)
    
    def __Calculate_AtmosLoss(self):
        """Calculates the Atmospheric Loss incurred by the signal using the ITU-R model.
        
        Returns: (double) AtmosLoss (dB)"""
        
        lon_GS, lat_GS = self.GS_Geo

        el = self.__Elevation_Angle
        f =  self.__Base_Frequency* u.GHz    # Link frequency
        p = 0.5 #percentage of time the value is exceeded
        Ag, Ac, Ar, As, A  =itur.atmospheric_attenuation_slant_path(lat_GS, lon_GS,
                                                                    f, el, p, D=5 * u.m,
                                                                    eta=0.65, return_contributions=True)

        self.__AtmosLoss = A.value

    def get_FSPL(self):
        return self.__FSPL
    
    def get_AtmosLoss(self):
        return self.__AtmosLoss

    