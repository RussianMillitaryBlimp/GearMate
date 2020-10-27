# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 17:24:06 2020

@author: Teak
"""

from __future__ import division
import numpy as np
import spur
# from matplotlib import pyplot as plt

#-----------------------------------------------------------------------------#
# Inputs

G_r = 3
D_c = 100e-3
PA = np.radians(20)
N_p = 20

standard_modulus_required = True # D_c will change to suit a standard module.
profile = "CUT" # CUT/MILLED, GROUND/SHAVED, CAST/FORGED, HOBBED/SHAPED
pinion_driven = 1 # If the pinion is driven, do calculations based on it.

T = 2 # Nm - Acting on PINION
om = 3000*(2*np.pi/60) # RPM -> rad/s
F =  5e-3 # in mm. Usually a factor of m.
v_p = 0.33
v_w = 0.33
E_p = 1e9
E_w = 1e9

#-----------------------------------------------------------------------------#
# Geometry

N_w, N_p = spur.geom_Gcalc(G_r, D_c, PA, N_p)
G_r = N_w/N_p

print("Adjusting gear ratio to %.2f"%G_r)

PCD = (D_c/(1+(1/G_r))*2, (D_c - D_c/(1+(1/G_r)))*2)

m = PCD[1]/N_p

p = np.pi*m

if standard_modulus_required == True:

    m = spur.choose_standard_m(m)

    PCD = (m*N_w, m*N_p)

    D_c = PCD[0]/2 + PCD[1]/2

    print("New centre distance of %.3f mm selected"%(D_c*1e3))

r_b = (np.cos(PA)*PCD[0]/2, np.cos(PA)*PCD[1]/2)

r_add = (m+PCD[0]/2, m+PCD[1]/2)

r_ded = (PCD[0]/2-1.25*m, PCD[1]/2-1.25*m)

t_root = spur.root_thickness(N_w, N_p, r_ded, r_b)

Lab = (r_add[0]**2 - r_b[0]**2)**0.5 + \
        (r_add[1]**2 - r_b[1]**2)**0.5 - \
        (PCD[0] + PCD[1])*np.sin(PA)/2

C_r = Lab/(p*np.cos(PA))

#-----------------------------------------------------------------------------#
# Check whether there is interference with current geometry

k = 1 # For fully formed teeth.

minN_p = np.ceil(2*k*(1 + (1 + 3*(np.sin(PA)**2))**0.5)/(3*(np.sin(PA)**2)))

if minN_p >= N_p:
    raise ValueError("Warning! Intersection detected between teeth.")

#-----------------------------------------------------------------------------#
# Forces

W_t = (2*T)/PCD[pinion_driven]
W_r = W_t*np.tan(PA)
W = (W_r**2 + W_t**2)**0.5

#-----------------------------------------------------------------------------#
# Stress

# Lewis bending stress

Kv = spur.Kv(om*PCD[pinion_driven]/2, profile)

Kf = spur.dolan_broghamer(min(t_root),
                          r_add[0] - r_ded[0],
                          r_ded[pinion_driven],
                          PCD[pinion_driven],
                          m,
                          PA)

Y = spur.lewis(t_root, r_add[0] - r_ded[0], m)

Sb = Kv*Kf*W_t/(F*m*Y[pinion_driven])

print("\nLewis bending stress is %.3f MPa"%(Sb*1e-6),
      "(facewidth of %.1f mm)"%(F*1e3))

# Hertzian contact stress

C_p = (1/(np.pi*((1 - (v_p**2))/E_p + (1 - (v_w)**2)/E_w)))**0.5

Sc = Kv*Kf*W_t*(2/(PCD[0]*np.sin(PA)) + 2/(PCD[1]*np.sin(PA)))/(F*np.cos(PA))

print("Herzian contact stress is %.3f MPa"%(-Sc*1e-6),
      "(facewidth of %.1f mm)"%(F*1e3))

#-----------------------------------------------------------------------------#
# Gear Life - This should be an extra section independant of the bending stress
#             calculations. The above should only give geometry and stresses,
#             and should steer clear of too much "engineers fudge".
#-----------------------------------------------------------------------------#
