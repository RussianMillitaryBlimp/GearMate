# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 17:24:06 2020

@author: Teak
"""

from __future__ import division
import numpy as np
import spur

#-----------------------------------------------------------------------------#
# Inputs

G_r = 3.33
D_c = 100e-3
PA = np.radians(20)
N_p = 20

standard_modulus_required = True # D_c will change to suit a standard module.
                                 # Note that the module will only be selected
                                 # between 0.5 < m < 50 (mm).

profile = "CUT" # CUT/MILLED, GROUND/SHAVED, CAST/FORGED, HOBBED/SHAPED

P = None  # kW
T = 2     # Nm
om = 3000 # RPM

F =  5e-3 # Facewidth in mm. Usually a factor of m.

v_p = 0.33
v_w = 0.33
E_p = 200e9
E_w = 200e9

#-----------------------------------------------------------------------------#
# Check whether there is interference with current geometry

minN_p = spur.minimum_involute_teeth(PA)

print("Minimum pinion tooth count to avoid interference is %d"%minN_p)

if minN_p >= N_p:
    raise ValueError("Warning! Intersection detected between teeth.")

#-----------------------------------------------------------------------------#
# Geometry

G_r, N = spur.adjust_gear_ratio(G_r, D_c, PA, N_p)

print("Adjusting gear ratio to %.2f for a whole tooth count."%G_r)

PCD = spur.define_PCD(G_r, D_c)

m, PCD, D_c, m_del = spur.choose_standard_m(N, PCD, standard_modulus_required)

print("Module of %.1f mm selected, m_delta = %.2f"%(m*1e3, m_del))
print("New centre distance of %.3f mm selected"%(D_c*1e3))

#-----------------------------------------------------------------------------#

pinion = spur.spurGear(PA, PCD[1], N[1], F)
pinion.E = E_p; pinion.v = v_p

wheel = spur.spurGear(PA, PCD[0], N[0], F)
wheel.E = E_w; wheel.v = v_w

C_r = pinion.contact_ratio(wheel)

#-----------------------------------------------------------------------------#

# Forces

W = pinion.force(T)

#-----------------------------------------------------------------------------#
# Stress

# Lewis bending stress

Sb = pinion.bending_stress(T, om)

print("\nLewis bending stress is %.3f MPa"%(Sb*1e-6),
      "(facewidth of %.1f mm)"%(pinion.facewidth*1e3))

# Hertzian contact stress

Sc = pinion.contact_stress(T, om, wheel)

print("Hertzian contact stress is %.3f MPa"%(-Sc*1e-6),
      "(facewidth of %.1f mm)"%(F*1e3))

#-----------------------------------------------------------------------------#
# Gear Life - This should be an extra section independant of the bending stress
#             calculations. The above should only give geometry and stresses,
#             and should steer clear of too much "engineers fudge".
#-----------------------------------------------------------------------------#
