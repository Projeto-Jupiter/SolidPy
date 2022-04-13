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
    def __init__(self, motor, propellant, initial_pressure=101325, empirical_data=None):
        self.motor = motor
        self.propellant = propellant
        self.initial_pressure = initial_pressure
        self.burn_area = self.motor.total_burn_area  # needs fixing
        self.set_parameters()
        self.time_span = 0.0, 2.5
        self.solve_burn()

        """Empirical known results (e.g. static-fire)"""
        self.empirical_time_steps, self.empirical_thrust = empirical_data
        self.empirical_chamber_pressure = self.empirical_evaluate_chamber_pressure()

    def set_parameters(self):
        self.parameters = (
            self.propellant.combustion_temperature,  # T_0
            self.propellant.products_constant,  # R
            self.propellant.density,  # rho_g
            self.propellant.specific_heat_ratio,  # k
            self.motor.nozzle_throat_area,  # A_t
        )

    def evaluate_nozzle_mass_flow(self, chamber_pressure):
        T_0, R, _, k, A_t = self.parameters
        return (
            chamber_pressure
            * A_t
            * np.sqrt(k / (R * T_0))
            * math.pow((2 / (k + 1)), ((k + 1) / (2 * (k - 1))))
        )

    def vector_field(self, time, state_variables):

        chamber_pressure, free_volume, regressed_length = state_variables
        T_0, R, rho_g, _, _ = self.parameters

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
            for (
                self.time,
                self.chamber_pressure,
                self.free_volume,
                self.regressed_length,
            ) in zip(
                self.solution.t,
                self.solution.y[0],
                self.solution.y[1],
                self.solution.y[2],
            ):
                self.evaluate_exit_pressure(self.chamber_pressure)
                self.evaluate_exit_velocity()
                self.evaluate_thrust(self.chamber_pressure)
                solution_writer.writerow(
                    [
                        self.time,
                        self.thrust,
                        self.chamber_pressure,
                        self.exit_pressure,
                        self.exit_velocity,
                        self.free_volume,
                        self.regressed_length,
                    ]
                )

    def evaluate_exit_mach(self):
        _, _, _, k, _ = self.parameters
        func = (
            lambda mach_number: math.pow((k + 1) / 2, -(k + 1) / (2 * (k - 1)))
            * math.pow((1 + (k - 1) / 2 * mach_number**2), (k + 1) / (2 * (k - 1)))
            / mach_number
            - self.motor.expansion_ratio
        )
        self.exit_mach = fsolve(func, 2)[0]
        return self.exit_mach

    def evaluate_exit_pressure(self, chamber_pressure):
        T_0, R, rho_g, k, A_t = self.parameters
        self.exit_pressure = chamber_pressure * math.pow(
            (1 + (k - 1) / 2 * self.evaluate_exit_mach() ** 2), -k / (k - 1)
        )
        return self.exit_pressure

    def evaluate_exit_velocity(self):
        T_0, R, rho_g, k, A_t = self.parameters
        self.exit_velocity = self.evaluate_exit_mach() * math.sqrt(k * R * T_0)
        return self.exit_velocity

    def evaluate_Cf(self):
        _, _, _, k, _ = self.parameters
        self.Cf = math.sqrt(
            (2 * k / (k - 1))
            * math.pow(2 / (k + 1), (k + 1) / (k - 1))
            * (1 - math.pow(self.exit_pressure / self.chamber_pressure, (k - 1) / k))
        )
        return self.Cf

    def evaluate_thrust(self, chamber_pressure):
        T_0, R, rho_g, k, A_t = self.parameters
        self.thrust = (
            self.evaluate_nozzle_mass_flow(chamber_pressure)
            * self.evaluate_exit_velocity()
            + (self.evaluate_exit_pressure(chamber_pressure) - self.initial_pressure)
            * self.motor.nozzle_exit_area
        )
        return self.thrust

    def evaluate_thrust2(self):
        self.evaluate_Cf()
        self.thrust = self.Cf * self.chamber_pressure * self.motor.nozzle_throat_area
        return self.thrust

    def empirical_evaluate_chamber_pressure(self):

        mach_number = self.evaluate_exit_mach()
        exit_velocity = self.evaluate_exit_velocity()
        T_0, R, rho_g, k, A_t = self.parameters

        chamber_pressure_list = []

        if self.empirical_thrust is not None:
            for thrust in self.empirical_thrust:
                func = lambda chamber_pressure: (
                    thrust
                    - self.evaluate_nozzle_mass_flow(chamber_pressure) * exit_velocity
                    + self.motor.nozzle_exit_area
                    * (
                        self.evaluate_exit_pressure(chamber_pressure)
                        - self.initial_pressure
                    )
                )

                current_chamber_pressure = fsolve(func, 25 * 10e5)
                chamber_pressure_list.append(current_chamber_pressure)

            return chamber_pressure_list

        return None

    def post_processing(self):
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

        figure(6, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            self.empirical_time_steps,
            self.empirical_thrust,
            "g",
            linewidth=0.75,
            label=r"$F^{emp}_T$",
        )
        legend(prop=FontProperties(size=16))
        title("Empirical Thrust as function of time")
        savefig("data/burn_simulation/graphs/empirical_thrust_pressure.png", dpi=200)

        figure(7, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            self.empirical_time_steps,
            self.empirical_chamber_pressure,
            "g",
            linewidth=0.75,
            label=r"$p^{emp}_c$",
        )
        legend(prop=FontProperties(size=16))
        title("Empirical Chamber Pressure as function of time")
        savefig("data/burn_simulation/graphs/empirical_chamber_pressure.png", dpi=200)

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

# Simulacao = Burn(Leviata, KNSB)

data_path = "data/static_fires/EmpuxoLeviata.csv"
ext_data = np.loadtxt(
    data_path,
    delimiter=",",
    unpack=True,
    skiprows=1,
    max_rows=2597,
    converters={0: lambda x: float((x[1:])[:-1]), 1: lambda y: float((y[1:])[:-1])},
)
#static fire data needs pre formatting

Simulacao = Burn(Leviata, KNSB, 101325, ext_data)

Simulacao.post_processing()
