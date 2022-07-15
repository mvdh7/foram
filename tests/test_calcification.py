# -*- coding: utf-8 -*-

#an ugly file o test things

#import foram
import numpy as np
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
step = 10
r = np.linspace(a,b,step) 

u = [np.ones(len(r))*C[0],
     np.ones(len(r))*C[1],
     np.ones(len(r))*C[2],
     np.ones(len(r))*C[3],
     np.ones(len(r))*C[4],
     np.ones(len(r))*C[5],
     np.ones(len(r))*C[6]]


sol = solve_bvp(odesys,bcs,r,u,t,s,a,Q,C)