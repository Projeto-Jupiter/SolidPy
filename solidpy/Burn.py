# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import csv
import math

import numpy as np
import scipy.constants as const

from scipy.optimize import fsolve
from scipy.integrate import solve_ivp, cumtrapz
from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, grid, legend, title, savefig, show

from Grain import Grain
from Propellant import Propellant
from Motor import Motor
from Environment import Environment
from Export import Export


class Burn:
    def __init__(self, grain, motor, propellant, environment=Environment()):
        self.motor = motor
        self.grain = grain
        self.propellant = propellant
        self.environment = environment

        self.gravity = environment.gravity
        self.environment_pressure = environment.atmospheric_pressure

        self.parameters = self.set_parameters()

    def set_parameters(self):
        parameters = (
            self.propellant.combustion_temperature,  # T_0
            self.propellant.products_constant,  # R
            self.propellant.density,  # rho_g
            self.propellant.specific_heat_ratio,  # k
            self.motor.nozzle_throat_area,  # A_t
        )
        return parameters

    def evaluate_nozzle_mass_flow(self, chamber_pressure):
        T_0, R, _, k, A_t = self.parameters
        return (
            chamber_pressure
            * A_t
            * np.sqrt(k / (R * T_0))
            * math.pow((2 / (k + 1)), ((k + 1) / (2 * (k - 1))))
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
        _, _, _, k, _ = self.parameters
        self.exit_pressure = chamber_pressure * math.pow(
            (1 + (k - 1) / 2 * self.evaluate_exit_mach() ** 2), -k / (k - 1)
        )
        return self.exit_pressure

    def evaluate_exit_temperature(self):
        T_0, _, _, k, _ = self.parameters
        self.exit_temperature = T_0 / (1 + (k - 1) / 2 * self.exit_mach**2)
        return self.exit_temperature

    def evaluate_exit_velocity(self):
        _, R, _, k, _ = self.parameters
        self.exit_velocity = self.evaluate_exit_mach() * math.sqrt(
            k * R * self.evaluate_exit_temperature()
        )
        return self.exit_velocity

    def evaluate_Cf(self, chamber_pressure):
        _, _, _, k, _ = self.parameters
        self.Cf = (
            math.sqrt(
                (2 * k**2 / (k - 1))
                * math.pow(2 / (k + 1), (k + 1) / (k - 1))
                * (
                    1
                    - math.pow(
                        self.evaluate_exit_pressure(chamber_pressure)
                        / chamber_pressure,
                        (k - 1) / k,
                    )
                )
            )
            + (
                (
                    self.evaluate_exit_pressure(chamber_pressure)
                    - self.environment_pressure
                )
                / chamber_pressure
            )
            * self.motor.expansion_ratio
        )
        return self.Cf

    def evaluate_thrust(self, chamber_pressure):
        self.thrust = (
            self.evaluate_Cf(chamber_pressure)
            * chamber_pressure
            * self.motor.nozzle_throat_area
        )
        return self.thrust

    def evaluate_total_impulse(self, thrust_list, time_list):
        total_impulse = cumtrapz(thrust_list, time_list)[-1]
        return total_impulse

    def evaluate_specific_impulse(self, thrust_list, time_list):
        specific_impulse = self.evaluate_total_impulse(thrust_list, time_list) / (
            self.propellant.density
            * self.grain.volume
            * self.motor.grain_number
            * self.gravity
        )
        return specific_impulse

    def evaluate_burn_rate(
        self, chamber_pressure, chamber_pressure_derivative, free_volume, burn_area
    ):
        T_0, R, rho_g, _, _ = self.parameters

        rho_0 = chamber_pressure / (R * T_0)  # product_gas_density
        nozzle_mass_flow = self.evaluate_nozzle_mass_flow(chamber_pressure)

        burn_rate = (
            free_volume / (R * T_0) * chamber_pressure_derivative + nozzle_mass_flow
        ) / (burn_area * (rho_g - rho_0))

        return burn_rate


class BurnSimulation(Burn):
    def __init__(self, grain, motor, propellant, environment=Environment()):
        Burn.__init__(self, grain, motor, propellant, environment)

        self.time_span = 0.0, 2.5
        self.number_steps = 1000

        self.solution = self.evaluate_solution()
        self.tail_off_solution = self.solve_tail_off_regime()

    """Solver required functions"""

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
            self.environment_pressure,
            self.motor.free_volume,
            regressed_length,
        ]

        def end_burn_propellant(time, state_variables):
            chamber_pressure, free_volume, regressed_length = state_variables
            if (
                (self.motor.chamber_volume - free_volume < 10e-6)
                or (31.92 / 2000 + regressed_length >= 71.92 / 2000)
                or (chamber_pressure == self.environment_pressure + 1)
            ):
                return 0
            return 1

        end_burn_propellant.terminal = True

        solution = solve_ivp(
            self.vector_field,
            self.time_span,
            state_variables,
            method="RK45",
            t_eval=np.linspace(*self.time_span, self.number_steps + 1),
            events=end_burn_propellant,
            # atol=absolute_error,
            # rtol=relative_error,
        )

        return solution

    def solve_tail_off_regime(self):

        T_0, R, _, _, A_t = self.parameters

        (
            self.initial_tail_off_time,
            self.initial_tail_off_chamber_pressure,
            self.initial_tail_off_free_volume,
        ) = (self.solution[0][-1], self.solution[1][-1], self.solution[2][-1])

        self.evaluate_tail_off_chamber_pressure = (
            lambda time: self.initial_tail_off_chamber_pressure
            * math.exp(
                -R
                * T_0
                * A_t
                / (self.initial_tail_off_free_volume * self.propellant.calc_cstar())
                * (time - self.initial_tail_off_time)
            )
        )

        self.tail_off_time = []
        self.tail_off_chamber_pressure = []

        for time in np.linspace(*self.time_span, self.number_steps + 1):
            tail_off_chamber_pressure = self.evaluate_tail_off_chamber_pressure(
                time + self.initial_tail_off_time
            )
            if tail_off_chamber_pressure / self.environment_pressure < 1.0001:
                break
            else:
                self.tail_off_time.append(time + self.initial_tail_off_time)
                self.tail_off_chamber_pressure.append(tail_off_chamber_pressure)

        self.tail_off_solution = [self.tail_off_time, self.tail_off_chamber_pressure]

        return self.tail_off_solution

    def evaluate_solution(self):
        burn_solution = self.solve_burn()

        thrust_list = []
        exit_velocity_list = []
        exit_pressure_list = []

        for chamber_pressure in burn_solution.y[0]:
            thrust_list.append(self.evaluate_thrust(chamber_pressure))
            exit_velocity_list.append(self.evaluate_exit_velocity())
            exit_pressure_list.append(self.evaluate_exit_pressure(chamber_pressure))

        burn_solution = [
            burn_solution.t,
            burn_solution.y[0],
            burn_solution.y[1],
            burn_solution.y[2],
            thrust_list,
            exit_pressure_list,
            exit_velocity_list,
        ]

        return burn_solution


class BurnExport(Export):
    def __init__(self, BurnSimulation):
        self.BurnSimulation = BurnSimulation
        self.solution = BurnSimulation.solution
        self.tail_off_solution = BurnSimulation.tail_off_solution

        self.post_processing()

    def post_processing(self):

        Export.raw_simulation_data_export(
            self.solution,
            "data/burn_simulation/burn_data.csv",
            [
                "Time",
                "Chamber Pressure",
                "Free Volume",
                "Regressed Length",
                "Thrust",
                "Exit Pressure",
                "Exit Velocity",
            ],
        )

        time, thrust = self.solution[0], self.solution[3]

        (
            self.max_chamber_pressure,
            self.end_free_volume,
            self.end_regressed_length,
            self.max_thrust,
            self.max_exit_pressure,
            self.max_exit_velocity,
        ) = Export.evaluate_max_variables_list(self.solution[0], self.solution[1:])

        self.total_impulse = self.BurnSimulation.evaluate_total_impulse(thrust, time)

        self.specific_impulse = self.BurnSimulation.evaluate_specific_impulse(
            thrust, time
        )

        return None

    def all_info(self):

        print("Max Thrust: ", self.max_thrust)
        print("Max Chamber pressure: ", self.max_chamber_pressure)
        print("Total impulse: ", self.total_impulse)
        print("Specific impulse: ", self.specific_impulse)

        return None

    def plotting(self):

        (
            time,
            chamber_pressure,
            free_volume,
            regressed_length,
            thrust,
            exit_pressure,
            _,
        ) = self.solution

        (tail_off_time, tail_off_chamber_pressure) = self.tail_off_solution

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

        figure(8, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            tail_off_time,
            tail_off_chamber_pressure,
            "m",
            linewidth=0.75,
            label=r"$p^{toff}_c$",
        )
        legend(prop=FontProperties(size=16))
        title("Tail Off Chamber Pressure as function of time")
        savefig("data/burn_simulation/graphs/tail_off_chamber_pressure.png", dpi=200)

        figure(9, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(time, chamber_pressure, "b", linewidth=0.75, label=r"$p_c$")
        legend(prop=FontProperties(size=16))
        plot(
            tail_off_time,
            tail_off_chamber_pressure,
            "m",
            linewidth=0.75,
            label=r"$p^{toff}_c$",
        )
        legend(prop=FontProperties(size=16))
        title("Total Burn Chamber Pressure as function of time")
        savefig("data/burn_simulation/graphs/total_burn_chamber_pressure.png", dpi=200)

        return None


if __name__ == "__main__":
    """Burn definitions"""
    Grao_Leviata = Grain(
        outer_radius=71.92 / 2000,
        initial_inner_radius=31.92 / 2000,
    )
    Leviata = Motor(
        Grao_Leviata,
        grain_number=4,
        chamber_inner_radius=77.92 / 2000,
        nozzle_throat_radius=17.5 / 2000,
        nozzle_exit_radius=44.44 / 2000,
        nozzle_angle=15 * np.pi / 180,
        chamber_length=600 / 1000,
    )
    KNSB = Propellant(
        specific_heat_ratio=1.1361,
        density=1700,
        products_molecular_mass=39.9e-3,
        combustion_temperature=1600,
        # burn_rate_a=5.13,
        # burn_rate_n=0.22,
        interpolation_list="data/burnrate/KNSB.csv",
    )
    Ambient = Environment(latitude=-0.38390456, altitude=750, ellipsoidal_model=True)

    """Static-fire data"""
    data_path = "data/static_fires/leviata_final_curve.csv"
    ext_data = np.loadtxt(
        data_path,
        delimiter=",",
        unpack=True,
        skiprows=1,
    )

    """Class instances"""
    Simulation = BurnSimulation(Grao_Leviata, Leviata, KNSB, Ambient)
    ExportPlot = BurnExport(Simulation)

    """Desired outputs"""
    ExportPlot.all_info()
    ExportPlot.plotting()
