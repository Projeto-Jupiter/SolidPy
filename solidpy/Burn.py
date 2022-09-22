# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = "MIT"
_license_ = ""

import math
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve
from scipy.integrate import solve_ivp, cumtrapz
from matplotlib.font_manager import FontProperties

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
        """Calculation of total nozzle mass flow.

        Source:
            https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html

        Args:
            chamber_pressure (float): current chamber pressure

        Returns:
            float: nozzle mass flow for the specified chamber pressure
        """
        T_0, R, _, k, A_t = self.parameters
        return (
            chamber_pressure
            * A_t
            * np.sqrt(k / (R * T_0))
            * math.pow((2 / (k + 1)), ((k + 1) / (2 * (k - 1))))
        )

    def evaluate_exit_mach(self):
        """Calculation of mach number at nozzle exit
        (ratio of flow speed to the local sound speed).

        Source:
        https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html

        Returns:
            float: mach number
        """
        _, _, _, k, _ = self.parameters
        func = (
            lambda mach_number: math.pow((k + 1) / 2, -(k + 1) / (2 * (k - 1)))
            * math.pow((1 + (k - 1) / 2 * mach_number**2), (k + 1) / (2 * (k - 1)))
            / mach_number
            - self.motor.expansion_ratio
        )
        self.exit_mach = fsolve(func, np.array(2))[0]
        return self.exit_mach

    def evaluate_exit_pressure(self, chamber_pressure):
        """Calculation of the pressure at nozzle exit .

        Source:
        https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html

        Args:
            chamber_pressure (float): current chamber pressure

        Returns:
            float: exit pressure for the specified chamber pressure
        """
        _, _, _, k, _ = self.parameters
        self.exit_pressure = chamber_pressure * math.pow(
            (1 + (k - 1) / 2 * self.evaluate_exit_mach() ** 2), -k / (k - 1)
        )
        return self.exit_pressure

    def evaluate_exit_temperature(self):
        """Calculation of fluid temperature at nozzle exit.

        Source:
        https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html

        Returns:
            float: exit temperature
        """
        T_0, _, _, k, _ = self.parameters
        self.exit_temperature = T_0 / (1 + (k - 1) / 2 * self.exit_mach**2)
        return self.exit_temperature

    def evaluate_exit_velocity(self):
        """Calculation of fluid velocity at nozzle exit.

        Source:
        https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html

        Returns:
            float: exit velocity
        """
        _, R, _, k, _ = self.parameters
        self.exit_velocity = self.evaluate_exit_mach() * math.sqrt(
            k * R * self.evaluate_exit_temperature()
        )
        return self.exit_velocity

    def evaluate_Cf(self, chamber_pressure):
        """Calculation of the engine's thrust coefficient.

        Source:
        Rogers, RasAero: The Solid Rocket Motor - Part 4 - Departures
        from Ideal Performance for Conical Nozzles and Bell Nozzles,
        page 28, eq.(9). High Power Rocketry.

        Args:
            chamber_pressure (float): current chamber pressure

        Returns:
            (float): the motor's thrust coefficient for
            a given chamber pressure
        """
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
        """Calculation of engine's thrust

        Args:
            chamber_pressure (float): current chamber

        Returns:
            float: motor's thrust for a given chamber pressure
        """
        self.thrust = (
            self.evaluate_Cf(chamber_pressure)
            * chamber_pressure
            * self.motor.nozzle_throat_area
        )
        return self.thrust

    def evaluate_total_impulse(self, thrust_list, time_list):
        """Numerical integration by trapezoids for total impulse
        approximation.

        Args:
            thrust_list (float list or float arrays): list of thrust values
            for each time step
            time_list (float list or float arrays): list of time steps

        Returns:
            float: the total impulse correspondent to the integral of
            the given values
        """
        total_impulse = cumtrapz(thrust_list, time_list)[-1]
        return total_impulse

    def evaluate_specific_impulse(self, thrust_list, time_list):
        """Calculation of motor's specific impulse.

        Args:
            thrust_list (float list or float arrays): list of thrust values
            time_list (float list or float arrays): list of time steps

        Returns:
            float: the specific impulse for the given values and propellant mass
        """
        specific_impulse = self.evaluate_total_impulse(thrust_list, time_list) / (
            self.propellant.density
            * self.grain.volume
            * self.motor.grain_number
            * self.environment.standard_gravity
        )
        return specific_impulse

    def evaluate_burn_rate(
        self, chamber_pressure, chamber_pressure_derivative, free_volume, burn_area
    ):
        """Calculation of propellant rate of regression, i.e. burn rate

        Args:
            chamber_pressure (float): current chamber pressure
            chamber_pressure_derivative (float): current chamber pressure derivative
            free_volume (float): current combustion chamber free volume
            burn_area (float): current total grain burn area accounting for regression

        Returns:
            float: current propellant burn rate
        """
        T_0, R, rho_g, _, _ = self.parameters

        rho_0 = chamber_pressure / (R * T_0)  # product_gas_density
        nozzle_mass_flow = self.evaluate_nozzle_mass_flow(chamber_pressure)

        burn_rate = (
            free_volume / (R * T_0) * chamber_pressure_derivative + nozzle_mass_flow
        ) / (burn_area * (rho_g - rho_0))

        return burn_rate


class BurnSimulation(Burn):
    def __init__(
        self,
        grain,
        motor,
        propellant,
        environment=Environment(),
        max_step_size=0.01,
        tail_off_evaluation=True,
    ):
        Burn.__init__(self, grain, motor, propellant, environment)
        self.max_step_size = max_step_size

        self.grain_burn_solution = self.evaluate_grain_burn_solution()
        self.tail_off_solution = (
            self.evaluate_tail_off_solution() if tail_off_evaluation else None
        )
        self.total_burn_solution = self.evaluate_complete_solution()

    """Solver required functions"""

    def vector_field(self, time, state_variables):
        """Generates the vector field of the corresponding simulations
        state variables (chamber pressure, free combustion chamber volume,
        length of grain regression), as required for solve_ivp differential
        equation solver. The grain regression length is measured with the
        zero at initial inner radius.

        Args:
            time (float): independent current time variable
            state_variables (list): simulation state variables
            to be solved

        Returns:
            list: vector field of the state variables
        """

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
        """Initial conditions setting and solver instantiation.

        Returns:
            object: solution object containing the solution
            for the differential equation state variables
        """
        regressed_length = 0
        state_variables = [
            self.environment_pressure,
            self.motor.free_volume,
            regressed_length,
        ]

        def end_burn_propellant(time, state_variables):
            """Establishment of solver terminal conditions. The simulation
            ends if the grain is totally regressed (burnt) emptying the
            combustion chamber.

            Args:
                time (float): independent current time variable
                state_variables
                (list): simulation state variables
                to be solved

            Returns:
                integer: boolean integer as termination parameter
            """
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
            method="DOP853",
            events=end_burn_propellant,
            max_step=0.01,
            atol=1e-8,
            rtol=1e-10,
        )

        return solution

    def solve_tail_off_regime(self):
        """Evaluates an analytical equation that describes the remaining
        chamber gases behavior after total grain burn.

        Returns:
            list: solution of the tail off regime, grouping time steps
            and chamber pressure
        """

        T_0, R, _, _, A_t = self.parameters

        # Set initial values at the end of grain burn simulation
        (
            self.initial_tail_off_time,
            self.initial_tail_off_chamber_pressure,
            self.initial_tail_off_free_volume,
        ) = (
            self.grain_burn_solution[0][-1],
            self.grain_burn_solution[1][-1],
            self.grain_burn_solution[2][-1],
        )

        # Analytical solution to the fluid behavior after grain burn
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

        # Keep the same time pacing for uniform union with burn solution
        time_steps = np.linspace(
            self.initial_tail_off_time,
            100.0,
            int((100.0 - self.initial_tail_off_time) / self.max_step_size),
        )

        tail_off_time = []
        tail_off_chamber_pressure = []
        tail_off_free_volume = []
        tail_off_regressed_length = []

        for time in time_steps:
            chamber_pressure = self.evaluate_tail_off_chamber_pressure(time)
            if chamber_pressure / self.environment_pressure > 1.0001:
                tail_off_time.append(time)
                tail_off_chamber_pressure.append(chamber_pressure)
                tail_off_free_volume.append(self.motor.chamber_volume)
                tail_off_regressed_length.append(
                    self.grain.outer_radius - self.grain.initial_inner_radius
                )
            else:
                break

        self.tail_off_solution = [
            tail_off_time,
            tail_off_chamber_pressure,
            tail_off_free_volume,
            tail_off_regressed_length,
        ]

        return self.tail_off_solution

    def process_solution(self, chamber_pressure_list):
        """Iteration through the solve_ivp solution in order to compute
        notable burn characteristics besides the state variables, such as thrust,
        exit pressure and exit velocity.

        Args:
            chamber_pressure_list (list): chamber pressure solution
            evaluated by solve_ivp

        Returns:
            tuple: thrust, exit pressure and exit velocity lists
        """
        thrust_list = []
        exit_velocity_list = []
        exit_pressure_list = []

        for chamber_pressure in chamber_pressure_list:
            thrust_list.append(self.evaluate_thrust(chamber_pressure))
            exit_velocity_list.append(self.evaluate_exit_velocity())
            exit_pressure_list.append(self.evaluate_exit_pressure(chamber_pressure))

        return thrust_list, exit_pressure_list, exit_velocity_list

    def evaluate_grain_burn_solution(self):
        """Adapts solve_ivp results to a simple matrix containing each
        burn characteristic.

        Returns:
            list: list containing the solution and added burn computations
        """
        grain_burn_solution = self.solve_burn()

        grain_burn_solution = [
            grain_burn_solution.t,
            grain_burn_solution.y[0],
            grain_burn_solution.y[1],
            grain_burn_solution.y[2],
            *self.process_solution(grain_burn_solution.y[0]),
        ]

        return grain_burn_solution

    def evaluate_tail_off_solution(self):
        """Adapts solve_ivp results to a simple matrix containing each
        burn characteristic.

        Returns:
            list: list containing the solution and added burn computations
        """
        tail_off_solution = self.solve_tail_off_regime()
        tail_off_solution = [
            *tail_off_solution,
            *self.process_solution(tail_off_solution[1]),
        ]

        return tail_off_solution

    def evaluate_complete_solution(self):
        """Groups both grain burn and tail-off solutions after processing.

        Returns:
            matrix: contains tail-off and grain burn regime
        """
        grain_solution = self.grain_burn_solution
        tail_off_solution = self.tail_off_solution
        total_burn_solution = []

        if tail_off_solution:
            for grain_parameter, tail_off_parameter in zip(
                grain_solution, tail_off_solution
            ):
                total_burn_solution.append(
                    np.append(grain_parameter, tail_off_parameter)
                )
        else:
            for grain_parameter in grain_solution:
                total_burn_solution.append(grain_parameter)

        return total_burn_solution


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
        ) = self.BurnSimulation.total_burn_solution

        self.burn_exporting()
        self.post_processing()

    def burn_exporting(self):
        """Method that calls Export class for solution exporting in a csv.

        Returns:
            None
        """
        try:
            Export.raw_simulation_data_export(
                self.BurnSimulation.total_burn_solution,
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
        except OSError as err:
            print("OS error: {0}".format(err))
        return None

    def post_processing(self):
        """Method for post solution values processing, allowing for final
        notable burn evaluations and notable solution points, such as extrema.

        Returns:
            None
        """
        (
            self.max_chamber_pressure,
            self.end_free_volume,
            self.end_regressed_length,
            self.max_thrust,
            self.max_exit_pressure,
            self.max_exit_velocity,
        ) = Export.evaluate_max_variables_list(
            self.BurnSimulation.total_burn_solution[0],
            self.BurnSimulation.total_burn_solution[1:],
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
        """Console logging of notable burn characteristics.

        Returns:
            None
        """
        print("Total Impulse: {:.2f} Ns".format(self.total_impulse))
        print("Max Thrust: {:.2f} N at {:.2f} s".format(*self.max_thrust))
        print("Mean Thrust: {:.2f} N".format(Export.evaluate_mean(self.thrust)))
        print(
            "Max Chamber Pressure: {:.2f} bar at {:.2f} s".format(
                self.max_chamber_pressure[0] / 1e5, self.max_chamber_pressure[1]
            )
        )
        print(
            "Mean Chamber Pressure: {:.2f} bar".format(
                Export.evaluate_mean(self.chamber_pressure) / 1e5
            )
        )
        print("Propellant mass: {:.2f} g".format(1000 * self.propellant_mass))
        print("Specific Impulse: {:.2f} s".format(self.specific_impulse))
        print("Burnout Time: {:.2f} s".format(self.time[-1]))

        return None

    def plotting(self):
        """Plot graphs of notable burn list values.

        Returns:
            None
        """
        plt.figure(1, figsize=(16, 9))
        plt.plot(self.time, self.thrust, color="b", linewidth=0.75, label=r"$F_T$")
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("thrust (N)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Thrust as function of time")
        plt.savefig("data/burn_simulation/graphs/thrust.png", dpi=200)

        plt.figure(2, figsize=(16, 9))
        plt.plot(
            self.time, self.chamber_pressure, color="b", linewidth=0.75, label=r"$p_c$"
        )
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("chamber pressure (pa)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Chamber Pressure as function of time")
        plt.savefig("data/burn_simulation/graphs/chamber_pressure.png", dpi=200)

        plt.figure(3, figsize=(16, 9))
        plt.plot(
            self.time, self.exit_pressure, color="b", linewidth=0.75, label=r"$p_e$"
        )
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("exit pressure (pa)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Exit Pressure as function of time")
        plt.savefig("data/burn_simulation/graphs/exit_pressure.png", dpi=200)

        plt.figure(4, figsize=(16, 9))
        plt.plot(
            self.time, self.free_volume, color="b", linewidth=0.75, label=r"$\forall_c$"
        )
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("free volume (m³)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Free Volume as function of time")
        plt.savefig("data/burn_simulation/graphs/free_volume.png", dpi=200)

        plt.figure(5, figsize=(16, 9))
        plt.plot(
            self.time,
            self.regressed_length,
            color="b",
            linewidth=0.75,
            label=r"$\ell_{regr}$",
        )
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("regressed length (m)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Regressed Grain Length as function of time")
        plt.savefig("data/burn_simulation/graphs/regressed_length.png", dpi=200)

        if self.BurnSimulation.tail_off_solution:
            plt.figure(6, figsize=(16, 9))
            plt.plot(
                self.BurnSimulation.tail_off_solution[0],
                self.BurnSimulation.tail_off_solution[1],
                color="b",
                linewidth=0.75,
                label=r"$p^{toff}_c$",
            )
            plt.grid(True)
            plt.xlabel("time (s)")
            plt.ylabel("tail of chamber pressure (pa)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Tail Off Chamber Pressure as function of time")
            plt.savefig(
                "data/burn_simulation/graphs/tail_off_chamber_pressure.png", dpi=200
            )

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
        # interpolation_list="data/burnrate/simulated/KNSB_Leviata_sim.csv",
    )

    Ambient = Environment(latitude=-0.38390456, altitude=627, ellipsoidal_model=True)

    """Class instances"""
    Simulation = BurnSimulation(
        Grao_Leviata, Leviata, KNSB, Ambient, tail_off_evaluation=True
    )
    ExportPlot = BurnExport(Simulation)

    """Desired outputs"""
    ExportPlot.all_info()
    ExportPlot.plotting()
