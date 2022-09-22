# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = "MIT"
_license_ = "x"

import math
import numpy as np
import scipy.interpolate as interpolate
import scipy.constants as const


class Propellant:
    def __init__(
        self,
        specific_heat_ratio,
        products_molecular_mass,
        combustion_temperature,
        density=None,
        **kwargs
    ):
        self.specific_heat_ratio = specific_heat_ratio
        self.products_molecular_mass = products_molecular_mass
        self.combustion_temperature = combustion_temperature
        self.density = density
        self.products_constant = const.R / products_molecular_mass
        self.cstar = self.evaluate_cstar()

        self.__dict__.update(kwargs)

    def evaluate_burn_rate(self, chamber_pressure):
        if "interpolation_list" in self.__dict__:
            pressure_list, burn_rate_list = np.loadtxt(
                self.interpolation_list, delimiter=",", unpack=True, skiprows=1
            )
            r_function = interpolate.interp1d(
                pressure_list, burn_rate_list, kind="cubic", fill_value="extrapolate"
            )
            r = r_function(chamber_pressure * 1e-6)  # r in mm/s; p in MPa
            return r / 1000

        elif "burn_rate_a" in self.__dict__ and "burn_rate_n" in self.__dict__:
            r = self.burn_rate_a * math.pow(chamber_pressure * 1e-6, self.burn_rate_n)
            return r / 1000

        raise TypeError(
            "Missing arguments. You must pass either an `interpolation_list` path or scalar ballistic "
            "coefficients `burn_rate_a` and `burn_rate_n` arguments to Propellant class "
        )

    def evaluate_cstar(self):
        k = self.specific_heat_ratio
        cstar = math.sqrt(
            (const.R * self.combustion_temperature)
            / (self.products_molecular_mass * k)
            * ((k + 1) / 2) ** ((k + 1) / (k - 1))
        )
        return cstar


if __name__ == "__main__":
    knsb = Propellant(
        1.1361,
        1700,
        39.86e-3,  # kg/mol
        1600,  # K
        interpolation_list=r"data/burnrate/KNSB.csv",
    )

    knsu = Propellant(
        1.1361, 1700, 39.86e-3, 1600, burn_rate_a=5.8, burn_rate_n=0.22  # kg/mol  # K
    )
