#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys 
sys.path.append("C:\\Users\\Bea\\Documents\\GitHub\\foram\\foram")
import numpy as np
from scipy.integrate import solve_bvp
import rates
import diffusion

# from https://scicomp.stackexchange.com/questions/36755/solve-a-system-of-coupled-differential-equations-in-python

def odesys(r,u): # maybe this needs to have the same input arguments as bcs, not sure
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
              (rates.k_p1(t,s)*co2 - rates.k_m1(t,s)*h*hco3 + rates.k_p4(t,s)*co2*oh - rates.k_m4(t,s)*hco3 + rates.k_p5(t,s)*h*co3 - rates.k_m5(t,s)*hco3),
            dco3, 
            -1/r * dco3 - 1/diffusion.Dc(t,"CO3",s)*
              (rates.k_m5(t,s)*hco3 - rates.k_p5(t,s)*h*co3),
            dh, 
            -1/r * dh - 1/diffusion.Dc(t,"H",s)*
              ((rates.k_m5(t,s) - rates.k_m1(t,s)*h)*hco3 + rates.k_p1(t,s)*co2 - rates.k_p5(t,s)*h*co3 + rates.k_p6(t,s) - rates.k_m6(t,s)*h*oh + rates.k_p7(t,s)*boh3 - rates.k_m7(t,s)*h*boh4),
            doh, 
            -1/r * oh - 1/diffusion.Dc(t,"OH",s)*
              (+ rates.k_m4(t,s)*hco3 - rates.k_p4(t,s)*co2*oh + rates.k_p6(t, s) - rates.k_m6(t,s)*h*oh)
            ,
            dboh3, 
            -1/r * dboh3 - 1/diffusion.Dc(t,"BOH3",s)*
              (-rates.k_p7(t,s)*boh3 + rates.k_m7(t,s)*h*boh4)
            ,
            dboh4, 
            -1/r * dboh4 - 1/diffusion.Dc(t,"BOH4",s)*
              (+ rates.k_p7(t,s)*boh3 - rates.k_m7(t,s)*h*boh4)
            ])


def bcs(ua,ub):
    """
    Evaluates residuals of the boundary conditions.
    For each of the 7 components of the carbonate system, we have a boundary condition
    for the first-order derivative of the concentration at r=a (fluxes; given by array Q) 
    and a boundary condition for the concentrations at r=b (given by array C)
    
    ua : values of concentrations and their first-order derivatives at r=a
    ub : values of concentrations and their first-order derivatives at r=b
    """
    #print(ua == u[:,1])
    return np.array([
        ub[0] - C[0],
        ua[1] - Q[0]/(4*np.pi*a**2 * diffusion.Dc_CO2(t)),
        ub[2] - C[1],
        ua[3] - Q[1]/(4*np.pi*a**2 * diffusion.Dc(t,"HCO3",s)),
        ub[4] - C[2],
        ua[5] - Q[2]/(4*np.pi*a**2 * diffusion.Dc(t,"CO3",s)),
        ub[6] - C[3],
        ua[7] - Q[3]/(4*np.pi*a**2 * diffusion.Dc(t,"H",s)),
        ub[8] - C[4],
        ua[9] - Q[4]/(4*np.pi*a**2 * diffusion.Dc(t,"OH",s)),
        ub[10] - C[5],
        ua[11] - Q[5]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH3",s)),
        ub[12] - C[6], 
        ua[13] - Q[6]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH4",s))
        ])



import PyCO2SYS as pyco2

#calculate the bulk conditions
#do pyco2sys
t = 24.5
s = 40.7

bulk = pyco2.sys(par1=2000, par2=8.25, par1_type= 2, par2_type=3, 
                 temperature=t, salinity=s)


C = [bulk["CO2"],
     bulk["HCO3"],
     bulk["CO3"],
     bulk["Hfree"],
     bulk["OH"],
     bulk["BOH3"],
     bulk["BOH4"]]

#calcification
Q = [0,0,3.25e-3,0,0,0,0]

a = 200
b = 10*a
r = np.linspace(a,b,11)

#make a matrix with initial guesses 
#that are not actually guesses but just the bulk concentrations
u = np.zeros((14, len(r)))
ion_names = np.array(["CO2", "HCO3", "CO3", "H", "OH", "BOH3", "BOH4"])
for i in range(14):
    for rad in range(len(r)):        
        if i%2 == 0:            
            u[i,rad] = C[int(i/2)]
        else:
             if i == 1: #for CO2
                 u[i, rad] = Q[0]/(4*np.pi*a**2 * diffusion.Dc_CO2(t))
             else: 
                 u[i, rad] = Q[int(np.floor(i/2))]/(4*np.pi*a**2 * diffusion.Dc(t,ion_names[int(np.floor(i/2))],s))
            
            
    
    

# u = [[C[0],C[0], C[0]],
#      [Q[0]/(4*np.pi*a**2 * diffusion.Dc_CO2(t)),Q[0]/(4*np.pi*a**2 * diffusion.Dc_CO2(t))],
#      [C[1],C[1]],
#      [Q[1]/(4*np.pi*a**2 * diffusion.Dc(t,"HCO3",s)),Q[1]/(4*np.pi*a**2 * diffusion.Dc(t,"HCO3",s))],
#      [C[2],92],
#      [Q[2]/(4*np.pi*a**2 * diffusion.Dc(t,"CO3",s)),Q[2]/(4*np.pi*a**2 * diffusion.Dc(t,"CO3",s))],
#      [C[3],C[3]],
#      [Q[3]/(4*np.pi*a**2 * diffusion.Dc(t,"H",s)),Q[3]/(4*np.pi*a**2 * diffusion.Dc(t,"H",s))],
#      [C[4],C[4]],
#      [Q[4]/(4*np.pi*a**2 * diffusion.Dc(t,"OH",s)),Q[4]/(4*np.pi*a**2 * diffusion.Dc(t,"OH",s))],
#      [C[5],C[5]],
#      [Q[5]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH3",s)),Q[5]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH3",s))],
#      [C[6],C[6]],
#      [Q[6]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH4",s)),Q[6]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH4",s))]]



sol = solve_bvp(odesys,bcs,r,u)



step = 2
#r_plot = np.arange(a,b+step,step)
co3_plot = sol.sol(r)[4]
import matplotlib.pyplot as plt
plt.figure()
plt.plot(r,co3_plot)
plt.scatter(r,u[4], marker="x")


# u = [[C[0],np.ones(len(r))*0,
#      np.ones(len(r))*C[1],np.ones(len(r))*0,
#      np.ones(len(r))*C[2],np.ones(len(r))*0,
#      np.ones(len(r))*C[3],np.ones(len(r))*0,
#      np.ones(len(r))*C[4],np.ones(len(r))*0,
#      np.ones(len(r))*C[5],np.ones(len(r))*0,
#      np.ones(len(r))*C[6],np.ones(len(r))*0],
#      [np.ones(len(r))*C[0],np.ones(len(r))*0,
#           np.ones(len(r))*C[1],np.ones(len(r))*0,
#           np.ones(len(r))*C[2],np.ones(len(r))*0,
#           np.ones(len(r))*C[3],np.ones(len(r))*0,
#           np.ones(len(r))*C[4],np.ones(len(r))*0,
#           np.ones(len(r))*C[5],np.ones(len(r))*0,
#           np.ones(len(r))*C[6],np.ones(len(r))*0]]
