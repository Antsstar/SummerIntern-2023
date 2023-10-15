import numpy as np
import math as m
import pandas as pd
import VectorClass as vc

class Terminal: #Terminal class for storing properties of a terminal.
    """ Define a User Terminal.

    Key Attributes:
    Geo - Stores the Longitude and Lattiude (currently assumes Elevation is zero but can easily be changed)
    Position - Stores the Current position of the terminal as a VectorClass object
    Pointing - Stores the pointing direction of the terminal, it should point towards the most 'optimal' satellite but can be normal to earth
    
    Key Methods:
    Properties(angle) - Performs a few functions to calculate the G/T, Gain, and EIRP of the Terminal for the given angle (degrees)
    """

    #Class attributes
    __Antenna_Noise_Temperature = 170 #K
    __Ambient_Temperature = 290 #K
    __Diplexer_Loss = 0.15 #dB
    __NF = 1.0 #dB
    __Peak_Gain = 33 #dB
    __Cosine_Roll_Off = 1.2


    def __init__(self, User, geo):
        self.__User = User
        self.__Geo = geo
        self.__Calc_Passive_Noise()
        self.__Calc_Active_Noise()
        

 ## Private funcitons   
    def __Calc_Passive_Noise(self):
        self.__Passive_Noise = self.__Ambient_Temperature * (10**(self.__Diplexer_Loss / 10) - 1)

    def __Calc_Active_Noise(self):
        self.__Active_Noise = self.__Ambient_Temperature * (10**(self.__NF / 10) - 1)

## Overrides
    def __str__(self):
        return (f"({self.__User}: {self.get_Geo()})")

## Setting functions
    def set_Antenna_Noise_Temperature(self, temp):
        self.__Antenna_Noise_Temperature = temp

    def set_Ambient_Temperature(self, temp):
        self.__Ambient_Temperature = temp
        self.__Calc_Active_Noise()
        self.__Calc_Passive_Noise()

    def set_Diplexer_Loss(self, temp):
        self.__Diplexer_Loss = temp
        self.__Calc_Passive_Noise()

    def set_NF(self, temp):
        self.__NF = temp
        self.__Calc_Active_Noise()

    def set_Position(self, array):
        """Takes an array of [x, y, z] cartesian coordiantes and defines a vector
        
        Parameters: (array) [x, y, z]
        """
        self.__Position = vc.Vector(*array, 'cartesian')
        self.__Pointing = vc.copy(self.__Position) #By default will point normal to Earth horizon

    def set_Pointing(self, array):
        """Takes an array of [x, y, z] cartesian coordiantes and defines a vector
        
        Parameters: (array) [x, y, z]
        """
        self.__Pointing = vc.Vector(*array, 'cartesian')

## Getting functions
    def get_Position(self):
        return self.__Position

    def get_Pointing(self):
        return self.__Pointing

    def get_Geo(self):
        return self.__Geo

# Calculate Noise Temperature function: Takes no inputs and outputs the total Noise temperature in dBK
    def __Calc_Noise_Temperature(self):
        Kelvin = self.__Antenna_Noise_Temperature + self.__Passive_Noise + self.__Active_Noise
        Noise_Temperature = 10 * m.log10(Kelvin)
        return Noise_Temperature

# Calculate Gain function: Takes one input, angle in degrees and outputs the terminal Gain in dBK
    def __Calc_Gain(self, angle):
        rad = m.radians(angle)
        Gain = self.__Peak_Gain + self.__Cosine_Roll_Off * 10 * m.log10(abs(m.cos(rad)))
        return Gain

    def Properties(self, angle):
        """Calculates the Terminal's Link Budget Properties.
        
        Paramters: Beam Angle in degrees
        
        Returns: (array) [G/T, G, EIRP]"""
        T = self.__Calc_Noise_Temperature()
        G = self.__Calc_Gain(angle)
        BUC = 9.03
        EIRP = BUC - self.__Diplexer_Loss + G
        return [G - T, G, EIRP]
