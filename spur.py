# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 17:24:06 2020

@author: Teak
"""

from __future__ import division
import numpy as np

def geom_Gcalc(G, D_c, PA, N_p):
    try:
        geom_Gcalc.call
    except:
        geom_Gcalc.call = False
    
    PCD_w = D_c/(1+(1/G))*2
    PCD_p = (D_c-PCD_w/2)*2
    
    p = np.pi*PCD_p/N_p
    N_w = int(np.pi*PCD_w/p)
    
    G = float(N_w/N_p)
    
    N_p = int(N_w/G)
    
    if geom_Gcalc.call == True:
        geom_Gcalc.call = False
        return N_w, N_p
    else:
        geom_Gcalc.call = True
        return geom_Gcalc(G, D_c, PA, N_p)


def choose_standard_m(m):
    standard = [0.5, 0.8, 1, 1.25, 1.5,
                2, 2.5, 3, 4, 5,
                6, 8, 10, 12, 16,
                20, 25, 32, 40, 50]
    m_new = min(standard, key=lambda x:abs(x-m*1e3))
    m_delta = m*1e3 - m_new
    
    print("Standard module of %.1f mm selected, "%m_new,
          "m_delta = %.2f"%m_delta)
    
    return m_new*1e-3


def inv(r, r_b):
    if r < r_b:
        PW = 0.
    else:
        PW = np.arccos(r_b/r)
    
    M = np.tan(PW) - PW
    
    return (r*np.cos(M), r*np.sin(M))


def inv_at_angle(f, phi, inverse = False):
    
    n = 1.
    
    if inverse == True:
        n = -1
    else:
        raise ValueError("Inverse input is incorect.")
    
    b = (f[0]*np.cos(phi)-n*f[1]*np.sin(phi),
         f[0]*np.sin(phi)+n*f[1]*np.cos(phi))
    
    return b


def root_thickness(N_w, N_p, r_ded, r_b):
    
    a = [None,None]; b = [None,None]; N = (N_w, N_p)
    
    for i in [0, 1]:
        a[i] = inv(r_ded[i], r_b[i])
        b[i] = inv_at_angle(a[i], 1.5*np.pi/N[i], True)
    
    l = (np.linalg.norm(np.array(a[0])-np.array(b[0])),
         np.linalg.norm(np.array(a[1])-np.array(b[1])))
    
    return l


def Kv(V, Profile = "CUT"):
    if Profile == "CUT":
        return (6.1 + V)/6.1
    elif Profile == "GROUND":
        return ((5.56 + V**0.5)/5.56)**0.5
    elif Profile == "CAST":
        return (3.05 + V)/3.05
    elif Profile == "HOBBED":
        return (3.56 + V**0.5)/3.56
    else:
        raise ValueError("Incorrect profile argument")


def lewis(t, l, m):
    return (t[0]**2 / (6 * m * l), t[1]**2 / (6 * m * l))

def dolan_broghamer(t, l, r_d, PCD, m, PA):
    r = (r_d - m)**2 / ((PCD/2) + r_d - m)  
    H = 0.340 - 0.4583662*PA  
    L = 0.316 - 0.4583662*PA
    M = 0.290 + 0.4583662*PA
    return H + ((t/r)**(L))*((t/l)**(M))
    
    
    