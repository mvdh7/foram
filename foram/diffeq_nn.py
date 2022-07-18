# -*- coding: utf-8 -*-

import elvet
#import sys 
#sys.path.append("C:\\Users\\Bea\\Documents\\GitHub\\foram\\foram")
import numpy as np
from scipy.integrate import solve_bvp
from foram import rates, diffusion

import PyCO2SYS as pyco2

#calculate the bulk conditions
#do pyco2sys
t = 24.5
s = 40.7

bulk = pyco2.sys(par1=2000, par2=8.25, par1_type= 2, par2_type=3, 
                 temperature=t, salinity=s)

#change units from umol/kg to mol/kg
C = np.array([bulk["CO2"],
     bulk["HCO3"],
     bulk["CO3"],
     bulk["Hfree"],
     bulk["OH"],
     bulk["BOH3"],
     bulk["BOH4"]])*1e-6

#calcification
#change units from nmol/h to mol/s
Q = np.array([0,0,3.25,0,0,0,0])*1e-9*3600

#change units from um to m
a = 200e-6
b = 100*a




def equations(r,c,dc):
    co2, hco3, co3, h, oh, boh3, boh4 = c[0], c[1], c[2], c[3], c[4], c[5], c[6]
    dco2, dhco3, dco3, dh, doh, dboh3, dboh4 = dc[0,0], dc[0,1], dc[0,2], dc[0,3], dc[0,4], dc[0,5], dc[0,6]
    return [
            -1/r * dco2 - 1/diffusion.Dc_CO2(t)*
              ((rates.k_m1(t,s)*h + rates.k_m4(t,s))*hco3 - (rates.k_p1(t,s) + rates.k_p4(t,s)*oh)*co2),
             
            -1/r * dhco3 - 1/diffusion.Dc(t,"HCO3",s)*
              (rates.k_p1(t,s)*co2 - rates.k_m1(t,s)*h*hco3 + rates.k_p4(t,s)*co2*oh - rates.k_m4(t,s)*hco3 + rates.k_p5(t,s)*h*co3 - rates.k_m5(t,s)*hco3),
             
            -1/r * dco3 - 1/diffusion.Dc(t,"CO3",s)*
              (rates.k_m5(t,s)*hco3 - rates.k_p5(t,s)*h*co3),
            
            -1/r * dh - 1/diffusion.Dc(t,"H",s)*
              ((rates.k_m5(t,s) - rates.k_m1(t,s)*h)*hco3 + rates.k_p1(t,s)*co2 - rates.k_p5(t,s)*h*co3 + rates.k_p6(t,s) - rates.k_m6(t,s)*h*oh + rates.k_p7(t,s)*boh3 - rates.k_m7(t,s)*h*boh4),
            
            -1/r * doh - 1/diffusion.Dc(t,"OH",s)*
              (+ rates.k_m4(t,s)*hco3 - rates.k_p4(t,s)*co2*oh + rates.k_p6(t, s) - rates.k_m6(t,s)*h*oh),
             
            -1/r * dboh3 - 1/diffusion.Dc(t,"BOH3",s)*
              (-rates.k_p7(t,s)*boh3 + rates.k_m7(t,s)*h*boh4),
             
            -1/r * dboh4 - 1/diffusion.Dc(t,"BOH4",s)*
              (+ rates.k_p7(t,s)*boh3 - rates.k_m7(t,s)*h*boh4)
            ]


bcs = [elvet.BC(b, lambda r, c, dc: c[0] - C[0]),
       elvet.BC(b, lambda r, c, dc: c[1] - C[1]),
       elvet.BC(b, lambda r, c, dc: c[2] - C[2]),
       elvet.BC(b, lambda r, c, dc: c[3] - C[3]),
       elvet.BC(b, lambda r, c, dc: c[4] - C[4]),
       elvet.BC(b, lambda r, c, dc: c[5] - C[5]),
       elvet.BC(b, lambda r, c, dc: c[6] - C[6]),
       elvet.BC(a, lambda r, c, dc: dc[0,0] - Q[0]/(4*np.pi*a**2 * diffusion.Dc_CO2(t))),
       elvet.BC(a, lambda r, c, dc: dc[0,1] - Q[1]/(4*np.pi*a**2 * diffusion.Dc(t,"HCO3",s))),
       elvet.BC(a, lambda r, c, dc: dc[0,2] - Q[2]/(4*np.pi*a**2 * diffusion.Dc(t,"CO3",s))),
       elvet.BC(a, lambda r, c, dc: dc[0,3] - Q[3]/(4*np.pi*a**2 * diffusion.Dc(t,"H",s))),
       elvet.BC(a, lambda r, c, dc: dc[0,4] - Q[4]/(4*np.pi*a**2 * diffusion.Dc(t,"OH",s))),
       elvet.BC(a, lambda r, c, dc: dc[0,5] - Q[5]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH3",s))),
       elvet.BC(a, lambda r, c, dc: dc[0,6] - Q[6]/(4*np.pi*a**2 * diffusion.Dc(t,"BOH4",s)))
       ]



domain = elvet.box((a , b, 100))

c = elvet.nn(1, 10, 7)

solver = elvet.solver(equations , bcs, domain, model=c , epochs =20000)

#%%
import elvet.plotting

elvet.plotting.plot_loss_density(solver)

