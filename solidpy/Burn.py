# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import math

import numpy as np

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
    def __init__(
        self, grain, motor, propellant, environment=Environment(), max_step_size=0.01
    ):
        Burn.__init__(self, grain, motor, propellant, environment)
        self.max_step_size = max_step_size

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
            if (self.motor.chamber_volume - free_volume < 1e-6) or (
                self.grain.inner_radius >= self.grain.outer_radius
            ):
                return 0
            return 1

        end_burn_propellant.terminal = True

        solution = solve_ivp(
            self.vector_field,
            (0.0, 100.0),
            state_variables,
            method="RK45",
            events=end_burn_propellant,
            max_step=0.01,
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

        time_steps = np.linspace(
            self.initial_tail_off_time,
            100.0,
            int((100.0 - self.initial_tail_off_time) / self.max_step_size),
        )

        self.tail_off_time = []
        self.tail_off_chamber_pressure = []

        for time in time_steps:
            tail_off_chamber_pressure = self.evaluate_tail_off_chamber_pressure(time)
            if tail_off_chamber_pressure / self.environment_pressure < 1.0001:
                break
            else:
                self.tail_off_time.append(time)
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
        (
            self.time,
            self.chamber_pressure,
            self.free_volume,
            self.regressed_length,
            self.thrust,
            self.exit_pressure,
            self.exit_velocity,
        ) = self.BurnSimulation.solution
        (
            self.tail_off_time,
            self.tail_off_chamber_pressure,
        ) = self.BurnSimulation.tail_off_solution

        self.post_processing()

    def post_processing(self):

        Export.raw_simulation_data_export(
            self.BurnSimulation.solution,
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

        (
            self.max_chamber_pressure,
            self.end_free_volume,
            self.end_regressed_length,
            self.max_thrust,
            self.max_exit_pressure,
            self.max_exit_velocity,
        ) = Export.evaluate_max_variables_list(
            self.BurnSimulation.solution[0], self.BurnSimulation.solution[1:]
        )

        self.total_impulse = self.BurnSimulation.evaluate_total_impulse(
            self.thrust, self.time
        )

        self.specific_impulse = self.BurnSimulation.evaluate_specific_impulse(
            self.thrust, self.time
        )

        self.propellant_mass = (
            self.BurnSimulation.motor.grain_number
            * self.BurnSimulation.motor.grain.volume
            * self.BurnSimulation.propellant.density
        )

        return None

    def all_info(self):

        print("Total Impulse: {:.2f} Ns".format(self.total_impulse))
        print("Max Thrust: {:.2f} N at {:.2f} s".format(*self.max_thrust))
        print("Mean Thrust: {:.2f} N".format(np.mean(self.thrust)))
        print(
            "Max Chamber Pressure: {:.2f} bar at {:.2f} s".format(
                self.max_chamber_pressure[0] / 1e5, self.max_chamber_pressure[1]
            )
        )
        print(
            "Mean Chamber Pressure: {:.2f} bar".format(
                np.mean(self.chamber_pressure) / 1e5
            )
        )
        print("Propellant mass: {:.2f} g".format(1000 * self.propellant_mass))
        print("Specific Impulse: {:.2f} s".format(self.specific_impulse))
        print("Burnout Time: {:.2f} s".format(self.time[-1]))

        return None

    def plotting(self):

        figure(1, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(self.time, self.thrust, "b", linewidth=0.75, label=r"$F_T$")
        legend(prop=FontProperties(size=16))
        title("Thrust as function of time")
        savefig("data/burn_simulation/graphs/thrust.png", dpi=200)

        figure(2, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(self.time, self.chamber_pressure, "b", linewidth=0.75, label=r"$p_c$")
        legend(prop=FontProperties(size=16))
        title("Chamber Pressure as function of time")
        savefig("data/burn_simulation/graphs/chamber_pressure.png", dpi=200)

        figure(3, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(self.time, self.exit_pressure, "b", linewidth=0.75, label=r"$p_e$")
        legend(prop=FontProperties(size=16))
        title("Exit Pressure as function of time")
        savefig("data/burn_simulation/graphs/exit_pressure.png", dpi=200)

        figure(4, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(self.time, self.free_volume, "b", linewidth=0.75, label=r"$\forall_c$")
        legend(prop=FontProperties(size=16))
        title("Free Volume as function of time")
        savefig("data/burn_simulation/graphs/free_volume.png", dpi=200)

        figure(5, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            self.time,
            self.regressed_length,
            "b",
            linewidth=0.75,
            label=r"$\ell_{regr}$",
        )
        legend(prop=FontProperties(size=16))
        title("Regressed Grain Length as function of time")
        savefig("data/burn_simulation/graphs/regressed_length.png", dpi=200)

        figure(8, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            self.tail_off_time,
            self.tail_off_chamber_pressure,
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
        plot(self.time, self.chamber_pressure, "b", linewidth=0.75, label=r"$p_c$")
        legend(prop=FontProperties(size=16))
        plot(
            self.tail_off_time,
            self.tail_off_chamber_pressure,
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
        interpolation_list="data/burnrate/KNSB3.csv",
    )
    Ambient = Environment(latitude=-0.38390456, altitude=627, ellipsoidal_model=True)

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
