# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

class Propellant:
    def __init__(
        self, 
        specific_heat_ratio, 
        density, 
        products_molecular_mass, 
        combustion_temperature,
        burn_rate_a,
        burn_rate_n,
        interpolation_list
    ):
        self.specific_heat_ratio = specific_heat_ratio
        self.density = density
        self.products_molecular_mass = products_molecular_mass
        # self.products_constant = R/products_molecular_mass
        self.combustion_temperature = combustion_temperature
        self.burn_rate_a = burn_rate_a
        self.burn_rate_n = burn_rate_n
        self.interpolation_list = interpolation_list
    
    def burn_rate(self, chamber_pressure):
        return #self.burn_rate_a * Math.exp(chamber_pressure, burn_rate_n)

    def mean_burn_rate(self):
        return

    def pressure_coeff(self):
        return
    
    def expontial_pressure_coeff(self):
        return

    def cstar(self):
        return 

