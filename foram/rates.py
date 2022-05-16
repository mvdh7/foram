import numpy as np
from . import constants


def k_p1(t, s=35):
    """Rate equation for CO2 + H2O => HCO3 + H (Table 1) in /s"""
    a1 = (680.5 - 4.72 * s) * 1e8  # /s
    e1 = 69.4  # kJ / mol
    return a1 * np.exp(-e1 / (constants.R * (t + 273.15)))

def k_m4(t, s=35):
    """Rate equation for CO2 + OH- -> HCO3- in /s"""
    return Kw * ((k_m1*k_p4) / k_p1)

def Kw (t, s=35):
    """Equilibrium constant for H2O -> H+ + OH-"""
    return np.exp(-13847.26/(t + 273.15) + 148.9652 - 23.6521 * np.log(t + 273.15) + (118.67/(t + 273.15) - 5.977 + 1.0495*np.log(t + 273.15))*s**0.5 - 0.01615 * s)
