import numpy as np
import math as m
import pandas as pd
import VectorClass as vc

class Satellite: #Satellite class for storing properties of a satellite.
    """ Define an Active Satellite

    Key Attributes:
    Geo - Stores Longitude, Lattitude and Elevation(km) of satellite
    Position - Stores the Current position of the satellite as a VectorClass object
    Pointing - Stores the pointing direction of the satellite, most of the time this will be Normal to the curvature of the 'Sphere'
    Name - Stores the name of the satellite
    EIRP - Stores transmitted power + Gain
    G/T - Calculated by manufacturer then provided. (changed)
    Downlink - Stores Downlink Frequency of satellite
    Uplink - Stores Uplink Frequency of satellite

    Key Methods:
    Properties(angle) - Performs a few functions to calculate the G/T, Gain, and EIRP of the Satellite for the given angle (degrees)
    """

    #Class Attributes
    __Peak_Gain = 33 #dB
    __Cosine_Roll_Off = 1.2


    def __init__(self, name, downlink=12.0, uplink=14.25, eirp=46.6, g_t=6):
        self.__Name = name
        self.__Downlink = downlink #GHz
        self.__Uplink = uplink #GHz
        self.__EIRP = eirp #dBW
        self.__G_T = g_t #dB/K


## Setting functions
    def set_Position(self, array, geo):
        """Takes an array of [x, y, z] cartesian coordiantes and defines a vector
        
        Parameters: (array) [x, y, z], (array) [Long, Lat, Elev]
        """
        self.__Position = vc.Vector(*array, "cartesian")
        self.__Geo = geo
        self.__Pointing =  -1 * self.__Position #By default points to Earth centre

    def set_Pointing(self, array):
        """Takes an array of [x, y, z] cartesian coordiantes and defines a vector
        
        Parameters: (array) [x, y, z]
        """
        self.__Pointing = vc.Vector(*array, "cartesian")


## Overrides
    def __str__(self):
        return (f"({self.__Name}, {self.get_Downlink()}, {self.get_Uplink()})")

## Getting functions
    def get_Uplink(self):
        return self.__Uplink

    def get_Downlink(self):
        return self.__Downlink

    def get_Position(self):
        return self.__Position

    def get_Pointing(self):
        return self.__Pointing    

    def get_Geo(self):
        return self.__Geo
    
    # Calculate Gain function: Takes one input, angle in degrees and outputs the terminal Gain in dBK
    def __Calc_Gain(self, angle):
        rad = m.radians(angle)
        Gain = self.__Peak_Gain + self.__Cosine_Roll_Off * 10 * m.log10(abs(m.cos(rad)))
        return Gain

    def Properties(self, angle):
        """Calculates the Satellite's Link Budget Properties.
        
        Paramters: Beam Angle in degrees
        
        Returns: (array) [G/T, G, EIRP]"""
        G_T = self.__G_T
        G = self.__Calc_Gain(angle)
        EIRP = self.__EIRP
        return [G_T, G, EIRP]
    

