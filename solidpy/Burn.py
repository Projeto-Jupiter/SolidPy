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
from scipy.integrate import solve_ivp
from scipy.integrate import cumtrapz
from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, grid, legend, title, savefig, show

from Grain import Grain
from Motor import Motor
from Propellant import Propellant



class Burn:
    def __init__(self, grain, motor, propellant, initial_pressure=101325, empirical_data=None):
        self.motor = motor
        self.grain = grain
        self.propellant = propellant
        self.initial_pressure = initial_pressure
        self.burn_area = self.motor.total_burn_area  # needs fixing
        self.set_parameters()
        self.time_span = 0.0, 2.5
        self.number_steps = 1000
        self.solve_burn()
        self.solve_tail_off_regime()
        self.raw_simulation_data_export()
        self.gravity = 9.81

        """Empirical known results (e.g. static-fire)"""
        self.empirical_time_steps, self.empirical_thrust = empirical_data
        self.empirical_chamber_pressure = self.empirical_evaluate_chamber_pressure()

    def set_parameters(self):
        self.parameters = (
            self.propellant.combustion_temperature,  # T_0
            self.propellant.products_constant,  # R
            self.evaluate_propellant_density(),  # rho_g
            self.propellant.specific_heat_ratio,  # k
            self.motor.nozzle_throat_area,  # A_t
        )
    
    def evaluate_propellant_density(self):
        if self.grain.density is not None:
            return self.grain.density
        return self.propellant.density

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

        def end_burn_propellant(time, state_variables):
            chamber_pressure, free_volume, regressed_length = state_variables
            if (
                (self.motor.chamber_volume - free_volume < 10e-6)
                or (31.92 / 2000 + regressed_length >= 71.92 / 2000)
                or (chamber_pressure == self.initial_pressure + 1)
            ):
                return 0
            return 1

        end_burn_propellant.terminal = True

        self.solution = solve_ivp(
            self.vector_field,
            self.time_span,
            state_variables,
            method="RK45",
            t_eval=np.linspace(*self.time_span, self.number_steps + 1),
            events=end_burn_propellant,
            # atol=absolute_error,
            # rtol=relative_error,
        )

        return None

    def solve_tail_off_regime(self):

        T_0, R, _, _, A_t = self.parameters

        (
            self.initial_tail_off_time,
            self.initial_tail_off_chamber_pressure,
            self.initial_tail_off_free_volume,
        ) = (self.solution.t[-1], self.solution.y[0][-1], self.solution.y[1][-1])

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
            if tail_off_chamber_pressure / self.initial_pressure < 1.0001:
                break
            else:
                self.tail_off_time.append(time + self.initial_tail_off_time)
                self.tail_off_chamber_pressure.append(tail_off_chamber_pressure)

        self.tail_off_solution = [self.tail_off_time, self.tail_off_chamber_pressure]

        return None

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

    def evaluate_Cf(self):
        _, _, _, k, _ = self.parameters
        self.Cf = math.sqrt(
            (2 * k / (k - 1))
            * math.pow(2 / (k + 1), (k + 1) / (k - 1))
            * (1 - math.pow(self.exit_pressure / self.chamber_pressure, (k - 1) / k))
        )
        return self.Cf

    def evaluate_thrust(self, chamber_pressure):
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

    # needs throughout revision, maybe better approach is integration scipy_cumtrapz
    def evaluate_specific_impulse(self, chamber_pressure_list, thrust_list):
        self.specific_impulse_list = []
        for chamber_pressure, thrust in zip(chamber_pressure_list, thrust_list):
            self.specific_impulse_list.append(
                thrust / (self.evaluate_nozzle_mass_flow(chamber_pressure))
            )
        specific_impulse = np.average(self.specific_impulse_list) / self.gravity
        return specific_impulse

    def empirical_evaluate_chamber_pressure(self):
        T_0, R, _, k, A_t = self.parameters
        chamber_pressure_list = []

        if self.empirical_thrust is not None:
            for thrust in self.empirical_thrust:
                current_chamber_pressure = (
                    thrust + self.motor.nozzle_exit_area * self.initial_pressure
                ) / (
                    A_t
                    * np.sqrt(k / (R * T_0))
                    * math.pow((2 / (k + 1)), ((k + 1) / (2 * (k - 1))))
                    * self.evaluate_exit_velocity()
                    + math.pow(
                        (1 + (k - 1) / 2 * self.evaluate_exit_mach() ** 2), -k / (k - 1)
                    )
                    * self.motor.nozzle_exit_area
                )

                chamber_pressure_list.append(current_chamber_pressure)

            return chamber_pressure_list

        return None

    def evaluate_max_list(self, list, time_list):
        max_index = np.argmax(list)
        max_value = list[max_index]
        time_of_max = time_list[max_index]
        return max_value, time_of_max

    def raw_simulation_data_export(self):

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

            return None

    def post_processing(self):

        self.simulation_data = np.loadtxt(
            "data/burn_simulation/burn_data.csv", delimiter=",", unpack=True, skiprows=1
        )
        (
            time,
            thrust,
            chamber_pressure,
            exit_pressure,
            exit_velocity,
            free_volume,
            regressed_length,
        ) = self.simulation_data

        def evaluate_max_variables_list():
            max_variable_list = []
            for data_variables in self.simulation_data[1:]:
                max_variable_list.append(self.evaluate_max_list(data_variables, time))
            max_variable_list.append(self.evaluate_max_list(self.empirical_thrust, self.empirical_time_steps))
            max_variable_list.append(
                self.evaluate_max_list(self.empirical_chamber_pressure, self.empirical_time_steps)
            )
            return max_variable_list

        (
            self.max_thrust,
            self.max_chamber_pressure,
            self.max_exit_pressure,
            self.max_exit_velocity,
            self.end_free_volume,
            self.end_regressed_length,
            self.max_empirical_thrust,
            self.max_empirical_chamber_pressure,
        ) = evaluate_max_variables_list()

        self.specific_impulse = self.evaluate_specific_impulse(chamber_pressure, thrust)
        self.empirical_specific_impulse = self.evaluate_specific_impulse(
            self.empirical_chamber_pressure, self.empirical_thrust
        )

        isp1 = cumtrapz(thrust, time)[-1] / (0.65 * 4 * 9.8)
        isp2 = cumtrapz(self.empirical_thrust, self.empirical_time_steps)[-1] / (0.65 * 4 * 9.8)
        print(isp1, isp2)

        return None

    def plotting(self):

        (
            time,
            thrust,
            chamber_pressure,
            exit_pressure,
            exit_velocity,
            free_volume,
            regressed_length,
        ) = self.simulation_data

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
        savefig("data/burn_simulation/graphs/empirical_thrust.png", dpi=200)

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

    def info():
        return None

    def all_info():
        return None


"""Rocket definitions"""
Grao_Leviata = Grain(outer_radius=71.92 / 2000, initial_inner_radius=31.92 / 2000, mass = 650 / 1000)
Leviata = Motor(
    Grao_Leviata,
    grain_number=4,
    chamber_inner_radius=77.92 / 2000,
    nozzle_throat_radius=17.5 / 2000,
    nozzle_exit_radius=44.44 / 2000,
    chamber_length=600 / 1000,
)
KNSB = Propellant(
    specific_heat_ratio=1.1361,
    density=1700,
    products_molecular_mass=39.9e-3,
    combustion_temperature=1600,
    burn_rate_a=5.13,
    burn_rate_n=0.22,
    # interpolation_list="data/burnrate/KNSB.csv",
)

"""Static-fire data"""
data_path = "data/static_fires/leviata_raw_data-1.csv"
ext_data = np.loadtxt(
    data_path,
    delimiter=";",
    unpack=True,
    skiprows=0,
)

Simulacao = Burn(Grao_Leviata, Leviata, KNSB, 101325, ext_data)

Simulacao.post_processing()
Simulacao.plotting()

print("Grain density: ", Simulacao.evaluate_propellant_density())
print("Max Thrust: ", Simulacao.max_thrust, Simulacao.max_empirical_thrust)
print("Max Chamber pressure: ", Simulacao.max_chamber_pressure, Simulacao.max_empirical_chamber_pressure)
print("Specific impulse: ", Simulacao.specific_impulse,Simulacao.empirical_specific_impulse)
