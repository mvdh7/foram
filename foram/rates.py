import numpy as np
from . import constants


def k_p1(t, s=35):
    """Rate equation for CO2 + H2O => HCO3 + H (Table 1) in /s"""
    a1 = (680.5 - 4.72 * s) * 1e8  # /s
    e1 = 69.4  # kJ / mol
    return a1 * np.exp(-e1 / (constants.R * (t + 273.15)))

# def k_p4(t, s=35):
       # """Rate equation for CO2 + HO- => HCO3- (Table 1) in kg /mol /s"""
a4 = 8718 # kJ/mol
e4 = 62.8 # kJ/mol
K4 = a4 * np.exp(-e4/(constants.R * (t + 273.15)))
