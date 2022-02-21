# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import csv
import numpy as np
import math
import scipy.constants as const

from Grain import Grain
from Motor import Motor
from Propellant import Propellant
from scipy.integrate import solve_ivp


class Burn:
    def __init__(self, motor, propellant, initial_pressure = 101325):
        self.motor = motor
        self.propellant = propellant
        self.initial_pressure = initial_pressure
        self.grain_area = motor.grain.burn_area
        self.set_parameters()
        self.time_span = 0.0, 2.0
        self.solve_burn()

    def set_parameters(self):
        self.parameters = (
            self.propellant.combustion_temperature,  # T_0
            self.propellant.products_constant,  # R
            self.propellant.density,  # rho_g
            self.propellant.specific_heat_ratio,  # k
            self.motor.nozzle_throat_area,  # A_t
        )

    def evaluate_nozzle_mass_flow(self, chamber_pressure):
        T_0, R, rho_g, k, A_t = self.parameters
        return (
            chamber_pressure
            * A_t
            * np.sqrt(k / (R * T_0))
            * math.pow((2 / (k + 1)), ((k + 1) / (2 * (k - 1))))
        )

    def vector_field(self, time, state_variables):

        chamber_pressure, free_volume, regressed_length = state_variables
        T_0, R, rho_g, k, A_t = self.parameters
        
        rho_0 = chamber_pressure / (R * T_0)  # product_gas_density
        nozzle_mass_flow = self.evaluate_nozzle_mass_flow(chamber_pressure)
        grain_area = self.motor.grain.evaluate_tubular_burn_area(regressed_length)
        burn_rate = self.propellant.burn_rate(chamber_pressure)
        print((grain_area * burn_rate * (rho_g - rho_0) - nozzle_mass_flow) * R * T_0 / free_volume)
        print(grain_area)
        print(rho_g)
        print(rho_0)
        print(const.R)
        print()

        vector_state = [
            (grain_area * burn_rate * (rho_g - rho_0) - nozzle_mass_flow) * R * T_0 / free_volume,
            grain_area * burn_rate,
            burn_rate
        ]

        return vector_state

    def solve_burn(self):
        regressed_length = 0
        state_variables = [self.initial_pressure, self.motor.free_volume, regressed_length]

        self.solution = solve_ivp(
            self.vector_field,
            self.time_span,
            state_variables,
            method="RK45",
            t_eval=np.linspace(0.0, 2.0, 100),
            # atol=absolute_error,
            # rtol=relative_error,
        )

        with open("data/burn_simulation/burn_data.csv", "w", newline='') as burn_data:
            solution_writer = csv.writer(burn_data)
            solution_writer.writerow(["Time", "Chamber Pressure", "Free Volume", "Regressed Length"])
            for solution_time, solution_chamber_pressure, solution_free_volume, solution_regressed_length in zip(
                self.solution.t, self.solution.y[0], self.solution.y[1], self.solution.y[2]
            ):
                if self.motor.chamber_volume - solution_free_volume < 10e-10:
                    break
                else:
                    solution_writer.writerow(
                        [solution_time, solution_chamber_pressure, solution_free_volume, solution_regressed_length]
                    )
                    
    def post_processing():
        return None

    def info():
        return None

    def all_info():
        return None

Grao_Leviata = Grain(outer_radius=71.92 / 2000, initial_inner_radius=31.92 / 2000)
Leviata = Motor(Grao_Leviata,grain_number=4,chamber_inner_radius=77.92 / 2000,nozzle_throat_radius=8.75 / 2000)
KNSB = Propellant(specific_heat_ratio=1.1361, density=1700, products_molecular_mass=39.86e-3, combustion_temperature=1600, burn_rate_a=5.8e-9, burn_rate_n=0.9, interpolation_list=None)
Simulacao = Burn(Leviata, KNSB, initial_pressure = 3000000)
print(Simulacao.solution.y[0])
