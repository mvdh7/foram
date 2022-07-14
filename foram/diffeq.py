#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import solve_bvp
from . import rates
from . import diffusion

# from https://scicomp.stackexchange.com/questions/36755/solve-a-system-of-coupled-differential-equations-in-python

def odesys(r,u,t,s=35): # maybe this needs to have the same input arguments as bcs, not sure
    """
    The system of ordinary differential equations (ODEs) describing the diffusion-reaction equations
    (Eqns. 23-29 of Wolf-Gladrow & Riesebell 1997)

    Parameters
    ----------
    r : radius
    u : array of concentrations and their first-order derivatives of 7 components of the carbonate system
    t : temperature
    s : salinity

    Returns
    -------
    array of first- and second-order derivatives of 7 components of the carbonate system

    """
    co2, dco2, hco3, dhco3, co3, dco3, h, dh, oh, doh, boh3, dboh3, boh4, dboh4 = u
    return np.array([dco2, 
            -1/r * dco2 - 1/diffusion.Dc_CO2(t)*
              ((rates.k_m1(t,s)*h + rates.k_m4(t,s))*hco3 - (rates.k_p1(t,s) + rates.k_p4(t,s)*oh)*co2),
            dhco3, 
            -1/r * dhco3 - 1/diffusion.Dc(t,"HCO3",s)*
              ...,
            dco3, 
            ,
            dh, 
            ,
            doh, 
            ,
            dboh3, 
            ,
            dboh4, ])


def bcs(ua,ub,t,s,a,Q,C):
    """
    Evaluates residuals of the boundary conditions.
    For each of the 7 components of the carbonate system, we have a boundary condition
    for the first-order derivative of the concentration at r=a (fluxes; given by array Q) 
    and a boundary condition for the concentrations at r=b (given by array C)
    
    ua : values of concentrations and their first-order derivatives at r=a
    ub : values of concentrations and their first-order derivatives at r=b
    """
    return np.array([
        ua[1] - Q[0]/(4*np.pi*a**2 * diffusion.Dc_CO2(t)),
        ua[3] - Q[1]/(4*np.pi*a**2 * diffusion.Dc(t,"HCO3",s)),
        ua[5] - Q[2]/(4*np.pi*a**2 * diffusion.Dc(t,"CO3",s)),
        ua[7] - Q[3]/(4*np.pi*a**2 * diffusion.Dc(t,"H",s)),
        ua[9] - Q[4]/(4*np.pi*a**2 * diffusion.Dc(t,"OH",s)),
        ua[11] - Q[5]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH3",s)),
        ua[13] - Q[6]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH4",s)),
        ub[0] - C[0],
        ub[2] - C[1],
        ub[4] - C[2],
        ub[6] - C[3],
        ub[8] - C[4],
        ub[10] - C[5],
        ub[12] - C[6]
        ])

a = 0
b = 1
step = 0.1
r = np.linspace(a,b,step)
u = np.array([])
Q = np.array([])
C = np.array([])

sol = solve_bvp(odesys,bcs,r,u,t,s,a,Q,C)
