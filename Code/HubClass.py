import numpy as np
import math as m
import pandas as pd
import VectorClass as vc

class Hub: #Hub class for storing properties of a terminal.
    """
    Key Attributes:
    Position - Stores the Current position of the terminal as a VectorClass object
    Pointing - Stores the pointing direction of the terminal, it should point towards the most 'optimal' satellite but can be normal to earth
    Peak_Gain (dBi) - Given by manufacturer
    Cosine_Roll_Off - Given by manufacturer
    Active Noise tempearture -
    Passive Noise temperature - Calculated ! Remember to store as dBK
    Diplexer loss - Given by manufacturer 0.15 dB
    Ambient Noise temperature - Temperature of surrondings usually around 290 K
    Antenna Noise temperature - Given by manufacturer 170 K for Kymeta!
    NF - (Noise figure) Given by manufacturer
    User - Stores the User information?

    """