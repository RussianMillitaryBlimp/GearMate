# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 17:24:06 2020

@author: Teak
"""

from __future__ import division
import numpy as np

#-----------------------------------------------------------------------------#
# Classes

class spurGear():
    def __init__(self, PA, PCD, N, F):
        self.pressureAngle = PA
        self.pitch = PCD
        self.teeth = N
        self.facewidth = F
        self.module = self.pitch/N
        self.addendum = self.pitch/2 + self.module
        self.dedendum = self.pitch/2 - 1.25*self.module
        self.base = np.cos(self.pressureAngle)*self.pitch/2
        self.root = spurGear.root_thickness(self)
        self.length = self.addendum - self.dedendum
        self.lewis = self.root**2 / (6 * self.module * self.length)
        self.Kf = spurGear.dolan_broghamer(self)
        
    def root_thickness(self):
        a = inv(self.dedendum, self.base)
        b = inv_at_angle(a, 1.5*np.pi/self.teeth, True)
        return np.linalg.norm(np.array(a)-np.array(b))
    
    def force(self, T):
        W_t = (2*T)/self.pitch
        W_r = W_t*np.tan(self.pressureAngle)
        W = (W_r**2 + W_t**2)**0.5
        self.forces = (W, W_r, W_t)
        return self.forces
    
    def dolan_broghamer(self):
        r = (self.dedendum - self.module)**2 / \
            ((self.pitch/2) + self.dedendum - self.module)
        H = 0.340 - 0.4583662*self.pressureAngle  
        L = 0.316 - 0.4583662*self.pressureAngle
        M = 0.290 + 0.4583662*self.pressureAngle
        return H + ((self.root/r)**(L))*((self.root/self.length)**(M))
    
    def action_length(self, gear):
        if type(gear) is not spurGear:
            raise TypeError("Argument is not a spurGear type.")
        elif round(self.pressureAngle, 3) != round(gear.pressureAngle, 3):
            raise ValueError("Pressure angle is different between gears.")
        else:        
            return (self.addendum**2 - self.base**2)**0.5 + \
                   (gear.addendum**2 - gear.base**2)**0.5 - \
                   (self.pitch + gear.pitch)*np.sin(self.pressureAngle)/2
    
    def contact_ratio(self, gear):
        return self.action_length(gear) / \
               (np.pi*self.module*np.cos(self.pressureAngle))
        
    def Kv(self, om, profile = "CUT"):
        V = self.pitch*om*np.pi/60
        if profile == "CUT":
            return (6.1 + V)/6.1
        elif profile == "GROUND":
            return ((5.56 + V**0.5)/5.56)**0.5
        elif profile == "CAST":
            return (3.05 + V)/3.05
        elif profile == "HOBBED":
            return (3.56 + V**0.5)/3.56
        else:
            raise ValueError("Incorrect profile argument")

    def bending_stress(self, T, om):
        self.Sb = self.Kf * self.Kv(om) * self.force(T)[2] / \
                  (self.facewidth * self.module * self.lewis)
        return self.Sb

    def contact_stress(self, T, om, gear):
        try:
            C_p = (1 / (np.pi*((1 - (self.v**2))/self.E + \
                               (1 - (gear.v**2))/gear.E)))**0.5
        except:
            raise AttributeError("Elastic gear properties undetected.")

        self.Sc = C_p * (self.Kf * self.Kv(om) * self.force(T)[2] * \
                  (2 / (self.pitch*np.sin(self.pressureAngle)) + \
                   2 / (gear.pitch*np.sin(gear.pressureAngle))) / \
                  (self.facewidth*np.cos(self.pressureAngle)))**0.5
        return self.Sc
    
    

#-----------------------------------------------------------------------------#
# Functions


def minimum_involute_teeth(PA):
    return np.ceil(2*(1 + (1 + 3*(np.sin(PA)**2))**0.5)/(3*(np.sin(PA)**2)))


def adjust_gear_ratio(G, D_c, PA, N_p):
    try:
        adjust_gear_ratio.call
    except:
        adjust_gear_ratio.call = False
    
    PCD = define_PCD(G, D_c)
    
    p = np.pi*PCD[1]/N_p
    N_w = int(np.pi*PCD[0]/p)
    
    G = float(N_w/N_p)
    
    N_p = int(N_w/G)
    
    if adjust_gear_ratio.call == True:
        adjust_gear_ratio.call = False
        return G, (N_w, N_p)
    else:
        adjust_gear_ratio.call = True
        return adjust_gear_ratio(G, D_c, PA, N_p)

def define_PCD(G, D_c):
    r = D_c/(1+(1/G))
    return (r*2, (D_c - r)*2)


def choose_standard_m(N, PCD, standard_modulus_required=True):
    
    standard = [0.5, 0.8, 1, 1.25, 1.5,
            2, 2.5, 3, 4, 5,
            6, 8, 10, 12, 16,
            20, 25, 32, 40, 50]
    
    m = PCD[1]/N[1]
    
    m_delta = 0
    
    if standard_modulus_required == True:
        
        m_new = min(standard, key=lambda x:abs(x-m*1e3))
        
        m_delta = m*1e3 - m_new
        
        m = m_new*1e-3
    
    PCD = (m*N[0], m*N[1])
    
    D_c = PCD[0]/2 + PCD[1]/2
    
    return m, PCD, D_c, m_delta


def inv(r, r_b):
    if r < r_b:
        PW = 0
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
