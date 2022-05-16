import numpy as np
from . import constants


def k_p1(t, s=35):
    """Rate equation for CO2 + H2O => HCO3 + H (Table 1) in /s"""
    a1 = (680.5 - 4.72 * s) * 1e8  # /s
    e1 = 69.4  # kJ / mol
    return a1 * np.exp(-e1 / (constants.R * (t + 273.15)))


def k_p5(t, s=35):
    """Rate equation for CO3 + H => HCO3 (Table 1) in kg /mol /s"""
    return 1e10  # kg /mol /s


def k_p6(t, s=35):
    """Rate equation for H2O => H + OH (Table 1) in mol /kg /s"""
    return 1.3e-3  # mol /kg /s

def k_m7(t, s=35):
    """Rate equation for BOH4 + H => BOH3 + H2O (Table 1) in kg /mol /s"""
    return 1e10  # kg /mol /s
