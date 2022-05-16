import numpy as np
from . import constants


def ionstr_DOE94(salinity):
    """Ionic strength following DOE94 (from PyCO2SYS)."""
    # === CO2SYS.m comments: =======
    # This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4.
    return 19.924 * salinity / (1000 - 1.005 * salinity)


def total_sulfate(s=35):
    """Total sulfate following DOE."""
    return 0.14 * s / (96.062 * 1.80655)


def kHSO4_FREE_D90a(t, s=35):
    """Bisulfate dissociation constant following D90a (from PyCO2SYS)."""
    # === CO2SYS.m comments: =======
    # Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
    # The goodness of fit is .021.
    # It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    # TYPO on p. 121: the constant e9 should be e8.
    # Output KS is on the free pH scale in mol/kg-sw.
    # This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
    TempK = t + 273.15
    Sal = s
    logTempK = np.log(TempK)
    IonS = ionstr_DOE94(Sal)
    lnKSO4 = (
        -4276.1 / TempK
        + 141.328
        - 23.093 * logTempK
        + (-13856 / TempK + 324.57 - 47.986 * logTempK) * np.sqrt(IonS)
        + (35474 / TempK - 771.54 + 114.723 * logTempK) * IonS
        + (-2698 / TempK) * np.sqrt(IonS) * IonS
        + (1776 / TempK) * IonS ** 2
    )
    return np.exp(lnKSO4) * (1 - 0.001005 * Sal)


def pH_free_to_total(t, s=35):
    """Free to Total pH scale conversion factor."""
    return 1.0 + total_sulfate(s=s) / kHSO4_FREE_D90a(t, s=s)


def pH_total_to_free(t, s=35):
    """Total to Free pH scale conversion factor."""
    return 1.0 / pH_free_to_total(t, s=s)


def K1(t, s=35):
    """Equilibrium constant for CO2 + H2O => HCO3 + H"""
    K1 = np.exp(
        (-2307.1266 / (t + 273.15))
        + 2.83655
        - 1.5529413 * np.log(t + 273.15)
        + ((-4.0484 / (t + 273.15)) - 0.20760841) * s ** (1 / 2)
        + 0.08468345 * s
        - 0.00654208 * s ** (3 / 2)
        + np.log(1 - 0.001005 * s)
    )  # mol /kg
    return K1 * pH_total_to_free(t, s=s)


def K2(t, s=35):
    """Equilibrium constant for HCO3- -> H+ + CO3--"""
    K2 = np.exp(
        -3351.610 / (t + 273.15)
        - 9.226508
        - 0.2005743 * np.log(t + 273.15)
        + (-23.9722 / (t + 273.15) - 0.106901773) * np.sqrt(s)
        + 0.1130822 * s
        - 0.00846934 * s ** 1.5
        + np.log(1 - 0.001005 * s)
    )
    return K2


def Kw(t, s=35):
    """Equilibrium constant for H2O -> H+ + OH-"""
    return np.exp(
        -13847.26 / (t + 273.15)
        + 148.9652
        - 23.6521 * np.log(t + 273.15)
        + (118.67 / (t + 273.15) - 5.977 + 1.0495 * np.log(t + 273.15)) * s ** 0.5
        - 0.01615 * s
    ) * pH_total_to_free(t, s=s)


def k_p1(t, s=35):
    """Rate equation for CO2 + H2O => HCO3 + H (Table 1) in /s"""
    a1 = (680.5 - 4.72 * s) * 1e8  # /s
    e1 = 69.4  # kJ / mol
    return a1 * np.exp(-e1 / (constants.R * (t + 273.15)))


def k_m1(t, s=35):
    """Rate equation for HCO3 + H => CO2 + H2O (Table 1) in kg /mol /s"""
    return k_p1(t, s=s) / K1(t, s=s)


def k_p4(t, s=35):
    """Rate equation for CO2 + HO- => HCO3- (Table 1) in kg /mol /s"""
    a4 = 8718  # kJ/mol
    e4 = 62.8  # kJ/mol
    return a4 * np.exp(-e4 / (constants.R * (t + 273.15)))


def k_m4(t, s=35):
    """Rate equation for CO2 + OH- -> HCO3- in /s"""
    return Kw(t, s=s) * k_m1(t, s=s) * k_p4(t, s=s) / k_p1(t, s=s)


def k_p5(t, s=35):
    """Rate equation for CO3 + H => HCO3 (Table 1) in kg /mol /s"""
    return 1e10  # kg /mol /s


def k_m5(t, s=35):
    """Rate equation for HCO3-  => CP3-- + H+ (Table 1) in /s"""
    return K2(t, s) * k_p5(t, s)


def k_p6(t, s=35):
    """Rate equation for H2O => H + OH (Table 1) in mol /kg /s"""
    return 1.3e-3  # mol /kg /s


def k_m6(t, s=35):
    """Rate equation for H + OH => H2O (Table 1) in kg /mol /s"""
    return k_p6(t, s=s) / Kw(t, s=s)


def k_m7(t, s=35):
    """Rate equation for BOH4 + H => BOH3 + H2O (Table 1) in kg /mol /s"""
    return 1e10  # kg /mol /s
