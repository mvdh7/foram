import numpy as np
from . import constants



def mu(t, s):
    
    """viscosity of freshwater, eq 12, reference: Siedler and Peters (1986)"""
    A = 1.1709
    B = 1.827e-3
    C = 89.93
    mu_20 = 1.002e-3 # N s /m2
    
    #not sure whether it is ln or log
    mu = np.exp((A*(20-t)- B*(t-20)**2)/(t + C))*mu_20
    
    if s > 0:
        """estimate viscosity of seawater, eq 14, based on data by Whitfield and Turner (1981)"""
        mu = mu/(0.9508 - 0.0007379*t)
       
    return mu
    

# TODO: add stuff for a
def Dc(t, a, s=35):
    """Diffusion, Stokes-Einstein relation), eq 11"""
    
    Dc = (constants.kB * (t+273.15)) / (6 * np.pi * mu(t, s) * a)
    
    return Dc
    
    
def Dc_CO2(t):
    """Diffusion coeff of CO2  by Jähne et al. (1987), eq 15"""
    a_co2 = 5019 * 1e-9 # m2 /s
    ea = 19.51    # kJ /mol
    
    return a_co2 * np.exp(-ea / (constants.R * (t+273.15))) 
