import numpy as np
from . import constants




def K1(t, s=35):
    """Equilibrium constant for CO2 + H2O => HCO3 + H"""
    K1 = np.exp((-2307.1266/(t+273.15)) + 2.83655 - 1.5529413 * np.log(t+273.15) 
                + ((-4.0484/(t+273.15)) - 0.20760841)* s**(1/2) + 0.08468345*s 
               - 0.00654208 * s**(3/2) + np.log(1-0.001005 * s)) # mol kg 
    return K1
  
  
def Kw (t, s=35):
    """Equilibrium constant for H2O -> H+ + OH-"""
    return np.exp(-13847.26/(t + 273.15) + 148.9652 - 23.6521 * np.log(t + 273.15) + (118.67/(t + 273.15) - 5.977 + 1.0495*np.log(t + 273.15))*s**0.5 - 0.01615 * s)


def k_p1(t, s=35):
    """Rate equation for CO2 + H2O => HCO3 + H (Table 1) in /s"""
    a1 = (680.5 - 4.72 * s) * 1e8  # /s
    e1 = 69.4  # kJ / mol
    return a1 * np.exp(-e1 / (constants.R * (t + 273.15)))


def k_m1(t, s=35):
    """Rate equation for HCO3 + H => CO2 + H2O (Table 1) in kg /mol /s"""
   
    return k_p1(t, s)/K1(t, s) 



def k_p4(t, s=35):
    """Rate equation for CO2 + HO- => HCO3- (Table 1) in kg /mol /s"""
    a4 = 8718 # kJ/mol
    e4 = 62.8 # kJ/mol
    return a4 * np.exp(-e4/(constants.R * (t + 273.15)))


def k_m4(t, s=35):
    """Rate equation for CO2 + OH- -> HCO3- in /s"""
    return Kw * ((k_m1*k_p4) / k_p1)




def k_p5(t, s=35):
    """Rate equation for CO3 + H => HCO3 (Table 1) in kg /mol /s"""
    return 1e10  # kg /mol /s


def k_p6(t, s=35):
    """Rate equation for H2O => H + OH (Table 1) in mol /kg /s"""
    return 1.3e-3  # mol /kg /s

def k_m7(t, s=35):
    """Rate equation for BOH4 + H => BOH3 + H2O (Table 1) in kg /mol /s"""
    return 1e10  # kg /mol /s
