# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.constants as const

from Grain import Grain
from Motor import Motor
from Propellant import Propellant
from scipy.integrate import solve_ivp
from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, grid, legend, title, savefig, show


class Burn:
    def __init__(self, motor, propellant, initial_pressure=101325):
        self.motor = motor
        self.propellant = propellant
        self.initial_pressure = initial_pressure
        self.burn_area = self.motor.total_burn_area  # needs fixing
        self.set_parameters()
        self.time_span = 0.0, 2.5
        self.solve_burn()
        self.post_processing()

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
        burn_area = (
            self.motor.grain_number
            * self.motor.grain.evaluate_tubular_burn_area(regressed_length)
        )
        burn_rate = self.propellant.evaluate_burn_rate(chamber_pressure)

        vector_state = [
            (burn_area * burn_rate * (rho_g - rho_0) - nozzle_mass_flow)
            * R
            * T_0
            / free_volume,
            burn_area * burn_rate,
            burn_rate,
        ]

        return vector_state

    def solve_burn(self):
        regressed_length = 0
        state_variables = [
            self.initial_pressure,
            self.motor.free_volume,
            regressed_length,
        ]

        def end_burn(time, state_variables):
            chamber_pressure, free_volume, regressed_length = state_variables
            if (
                (self.motor.chamber_volume - free_volume < 10e-10)
                or (31.92 / 2000 + regressed_length >= 71.92 / 2000)
                or (chamber_pressure == self.initial_pressure + 1)
            ):
                return 0
            return 1

        end_burn.terminal = True

        self.solution = solve_ivp(
            self.vector_field,
            self.time_span,
            state_variables,
            method="RK45",
            t_eval=np.linspace(0.0, 2.5, 1001),
            events=end_burn,
            # atol=absolute_error,
            # rtol=relative_error,
        )

        with open("data/burn_simulation/burn_data.csv", "w", newline="") as burn_data:
            solution_writer = csv.writer(burn_data)
            solution_writer.writerow(
                [
                    "Time",
                    "Thrust",
                    "Chamber Pressure",
                    "Exit Pressure",
                    "Exit Velocity",
                    "Free Volume",
                    "Regressed Length",
                ]
            )
            for time, chamber_pressure, free_volume, regressed_length in zip(
                self.solution.t,
                self.solution.y[0],
                self.solution.y[1],
                self.solution.y[2],
            ):
                exit_pressure = self.evaluate_exit_pressure(chamber_pressure)
                exit_velocity = self.evaluate_exit_velocity()
                thrust = self.evaluate_thrust(chamber_pressure)
                solution_writer.writerow(
                    [
                        time,
                        thrust,
                        chamber_pressure,
                        exit_pressure,
                        exit_velocity,
                        free_volume,
                        regressed_length,
                    ]
                )

    def evaluate_exit_mach(self):
        T_0, R, rho_g, k, A_t = self.parameters
        func = (
            lambda mach_number: math.pow((k + 1) / 2, -(k + 1) / (2 * (k - 1)))
            * math.pow((1 + (k - 1) / 2 * mach_number**2), (k + 1) / (2 * (k - 1)))
            / mach_number
            - self.motor.expansion_ratio
        )
        mach = fsolve(func, 2)
        return mach

    def evaluate_exit_pressure(self, chamber_pressure):
        T_0, R, rho_g, k, A_t = self.parameters
        exit_pressure = chamber_pressure * math.pow(
            (1 + (k - 1) / 2 * self.evaluate_exit_mach()[0] ** 2), -k / (k - 1)
        )
        return exit_pressure

    def evaluate_exit_velocity(self):
        T_0, R, rho_g, k, A_t = self.parameters
        exit_velocity = self.evaluate_exit_mach()[0] * math.sqrt(k * R * T_0)
        return exit_velocity

    def evaluate_thrust(self, chamber_pressure):
        T_0, R, rho_g, k, A_t = self.parameters
        thrust = (
            self.evaluate_nozzle_mass_flow(chamber_pressure)
            * self.evaluate_exit_velocity()
            + (self.evaluate_exit_pressure(chamber_pressure) - self.initial_pressure)
            * self.motor.nozzle_exit_area
        )
        return thrust

    def post_processing(self):
        print(self.evaluate_exit_mach())
        print(self.evaluate_exit_pressure(4073878.314629354))
        print(self.evaluate_exit_pressure(4073878.314629354) / self.initial_pressure)
        print(self.evaluate_thrust(4073878.314629354))

        (
            time,
            thrust,
            chamber_pressure,
            exit_pressure,
            exit_velocity,
            free_volume,
            regressed_length,
        ) = np.loadtxt(
            "data/burn_simulation/burn_data.csv", delimiter=",", unpack=True, skiprows=1
        )

        figure(1, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(time, thrust, "b", linewidth=0.75, label=r"$F_T$")
        legend(prop=FontProperties(size=16))
        title("Thrust as function of time")
        savefig("data/burn_simulation/graphs/thrust.png", dpi=200)

        figure(2, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(time, chamber_pressure, "b", linewidth=0.75, label=r"$p_c$")
        legend(prop=FontProperties(size=16))
        title("Chamber Pressure as function of time")
        savefig("data/burn_simulation/graphs/chamber_pressure.png", dpi=200)

        figure(3, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(time, exit_pressure, "b", linewidth=0.75, label=r"$p_e$")
        legend(prop=FontProperties(size=16))
        title("Exit Pressure as function of time")
        savefig("data/burn_simulation/graphs/exit_pressure.png", dpi=200)

        figure(4, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(time, free_volume, "b", linewidth=0.75, label=r"$\forall_c$")
        legend(prop=FontProperties(size=16))
        title("Free Volume as function of time")
        savefig("data/burn_simulation/graphs/free_volume.png", dpi=200)

        figure(5, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(time, regressed_length, "b", linewidth=0.75, label=r"$\ell_{regr}$")
        legend(prop=FontProperties(size=16))
        title("Regressed Grain Length as function of time")
        savefig("data/burn_simulation/graphs/regressed_length.png", dpi=200)

        return None

    def info():
        return None

    def all_info():
        return None


Grao_Leviata = Grain(outer_radius=71.92 / 2000, initial_inner_radius=31.92 / 2000)
Leviata = Motor(
    Grao_Leviata,
    grain_number=4,
    chamber_inner_radius=77.92 / 2000,
    nozzle_throat_radius=17.5 / 2000,
    nozzle_exit_radius=44.44 / 2000,
)
KNSB = Propellant(
    specific_heat_ratio=1.1361,
    density=1700,
    products_molecular_mass=39.9e-3,
    combustion_temperature=1600,
    interpolation_list="data/burnrate/KNSB.csv",
)
Simulacao = Burn(Leviata, KNSB)
