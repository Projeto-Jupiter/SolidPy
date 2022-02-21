# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = "x"
_license_ = "x"

import math
import scipy.constants as const


class Propellant:
    def __init__(
        self,
        specific_heat_ratio,
        density,
        products_molecular_mass,
        combustion_temperature,
        burn_rate_a,
        burn_rate_n,
        interpolation_list,
    ):
        self.specific_heat_ratio = specific_heat_ratio
        self.density = density
        self.products_molecular_mass = products_molecular_mass
        self.products_constant = const.R / products_molecular_mass
        self.combustion_temperature = combustion_temperature
        self.burn_rate_a = burn_rate_a
        self.burn_rate_n = burn_rate_n
        self.interpolation_list = interpolation_list
        self.calc_cstar()

        self.mean_burn_rate_index = 0
        self.mean_burn_rate_value = 0

    def burn_rate(self, chamber_pressure):
        r = self.burn_rate_a * math.pow(chamber_pressure, self.burn_rate_n)
        self.mean_burn_rate(r)
        return r

    def calc_cstar(self):
        k = self.specific_heat_ratio
        self.cstar = math.sqrt(
            (const.R * self.combustion_temperature)
            / (self.products_molecular_mass * k)
            * ((k + 1) / 2) ** ((k + 1) / (k - 1))
        )
        return self.cstar

    # Test - Média Móvel
    def mean_burn_rate(self, current_burn_rate):
        self.mean_burn_rate_value *= self.mean_burn_rate_index
        self.mean_burn_rate_value += current_burn_rate
        self.mean_burn_rate_index += 1
        self.mean_burn_rate_value /= self.mean_burn_rate_index

    def pressure_coeff(self):
        return

    def exponential_pressure_coeff(self):
        return


knsb = Propellant(
    1.1361,
    1700,
    39.86e-3, #kg/mol
    1600, #K
    5.8e-9,
    0.9,
    None
)

