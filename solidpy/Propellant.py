# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = "x"
_license_ = "x"

import math
import csv
import scipy.interpolate as interpolate
import scipy.constants as const


class Propellant:
    def __init__(
        self,
        specific_heat_ratio,
        products_molecular_mass,
        combustion_temperature,
        density = None,
        **kwargs
    ):
        self.specific_heat_ratio = specific_heat_ratio
        self.density = density
        self.products_molecular_mass = products_molecular_mass
        self.products_constant = const.R / products_molecular_mass
        self.combustion_temperature = combustion_temperature

        self.__dict__.update(kwargs)
        self.calc_cstar()

        self.mean_burn_rate_index = 0
        self.mean_burn_rate_value = 0

    def evaluate_burn_rate(self, chamber_pressure):
        if "interpolation_list" in self.__dict__:
            with open(self.interpolation_list, "r") as interpolation_data:  # open file
                reader = csv.reader(interpolation_data)  # read csv
                next(reader)  # skip header
                # initialize and store data
                burn_rate_list = []
                pressure_list = []
                for line in reader:
                    burn_rate_list.append(float(line[1]))
                    pressure_list.append(float(line[0]))
            r_function = interpolate.interp1d(
                pressure_list, burn_rate_list, kind="cubic"
            )
            if chamber_pressure < 0:
                return 0
            r = r_function(chamber_pressure * 1e-6)  # r in mm/s; p in MPa
            return r / 1000

        elif "burn_rate_a" in self.__dict__ and "burn_rate_n" in self.__dict__:
            r = self.burn_rate_a * math.pow(chamber_pressure * 1e-6, self.burn_rate_n)
            return r / 1000

        else:
            raise TypeError(
                "Missing arguments. You must pass either an `interpolation_list` path or scalar ballistic coefficients `burn_rate_a` and `burn_rate_n` arguments to Propellant class"
            )

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


# knsb = Propellant(
#     1.1361,
#     1700,
#     39.86e-3, #kg/mol
#     1600, #K
#     interpolation_list=r'C:\Users\ansys\Desktop\SolidPy\data\burnrate\KNSB.csv'
# )

# knsu = Propellant(
#     1.1361,
#     1700,
#     39.86e-3, #kg/mol
#     1600, #K
#     burn_rate_a=5.8,
#     burn_rate_n=0.22
# )
