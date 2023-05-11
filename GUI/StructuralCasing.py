# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo Dos Santos De Souza"
_copyright_ = "MIT"
_license_ = ""

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve
from scipy.integrate import solve_ivp, cumtrapz
from matplotlib.font_manager import FontProperties


class Material:
    def __init__(
        self,
        ultimate_tensile_strength,
        shear_strength,
        yield_tensile_strength,
        bearing_strength,
    ):
        self.UTS = ultimate_tensile_strength
        self.shear_strenght = shear_strength
        self.YTS = yield_tensile_strength
        self.BYS = bearing_strength

class Geometry:
    def __init__(
        self,
        casing_inside_diameter,
        minor_bolt_diameter,
        bolts_number,
        wall_thickness,
        edge_distance,
        major_bolt_diameter
    ):
        self.casing_inside_diameter = casing_inside_diameter
        self.minor_bolt_diameter = minor_bolt_diameter
        self.bolts_number = bolts_number
        self.wall_thickness = wall_thickness
        self.edge_distance = edge_distance
        self.casing_outer_diameter = casing_inside_diameter + 2*wall_thickness
        self.major_bolt_diameter = major_bolt_diameter

    def evaluate_E_min(self):
        min_edge_distance = self.edge_distance - self.major_bolt_diameter/2

        return min_edge_distance

class StructuralAnalysis:
    def __init__(self,MEOP,material,geometry):
        self.MEOP = MEOP
        self.material = material
        self.geometry = geometry

    def bolt_shear(self):
        stress = (self.geometry.casing_inside_diameter**2)*self.MEOP/(self.geometry.bolts_number*self.geometry.minor_bolt_diameter**2)
        safety_factor = 0.75*self.material.UTS/stress
        
        return safety_factor

    def bolt_tear_out(self):
        f_bolt = (np.pi/4)*(self.geometry.casing_inside_diameter**2)*self.MEOP/(self.geometry.bolts_number)
        stress = f_bolt/(self.geometry.evaluate_E_min()*2*self.geometry.wall_thickness)
        safety_factor = self.material.shear_strenght/stress

        return safety_factor
    
    def casing_tensile(self):
        stress = (np.pi/4)*(self.geometry.casing_inside_diameter**2)*self.MEOP/(((self.geometry.casing_outer_diameter-self.geometry.wall_thickness)*np.pi-self.geometry.bolts_number*self.geometry.major_bolt_diameter)*self.geometry.wall_thickness)
        safety_factor = self.material.YTS / stress

        return safety_factor

    def bearing(self):
        f_bolt = (np.pi/4)*(self.geometry.casing_inside_diameter**2)*self.MEOP/(self.geometry.bolts_number)
        stress = f_bolt/(self.geometry.major_bolt_diameter*self.geometry.wall_thickness)
        safety_factor = self.material.BYS/stress

        return(safety_factor)

    def casing_circumferential_von_misses(self):
        sigma_1 = self.MEOP*self.geometry.casing_inside_diameter/(2*self.geometry.wall_thickness)
        sigma_2 = self.MEOP*self.geometry.casing_inside_diameter/(4*self.geometry.wall_thickness)
        equivalent_stress = np.sqrt(((sigma_1-sigma_2)**2+(sigma_1)**2+(sigma_2)**2)/2)

        safety_factor = self.material.YTS / equivalent_stress

        return safety_factor

