import numpy as np
import PyCO2SYS as pyco2
from . import constants


def k_p1(t, s=35):
    """Rate equation for CO2 + H2O => HCO3 + H (Table 1) in /s"""
    a1 = (680.5 - 4.72 * s) * 1e8  # /s
    e1 = 69.4  # kJ / mol
    return a1 * np.exp(-e1 / (constants.R * (t + 273.15)))

def k_m1(t, s=35):
    """Rate equation for HCO3 + H => CO2 + H2O (Table 1) in kg /mol /s"""
   
    return k_p1(t, s)/K1(t, s) 


def K1(t, s=35):
    """Equilibrium constant for CO2 + H2O => HCO3 + H"""
    
    K1 = np.exp((-2307.1266/(t+273.15)) + 2.83655 - 1.5529413 * np.log(t+273.15) 
                + ((-4.0484/(t+273.15)) - 0.20760841)* s**(1/2) + 0.08468345*s 
               - 0.00654208 * s**(3/2) + np.log(1-0.001005 * s)) # mol kg 
    
    #K1 = pyco2.k_carbonic_1(temperature=t, salinity=2, opt_k_carbonic=1)
    
    return K1

