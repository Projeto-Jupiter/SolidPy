# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.integrate import solve_ivp
from matplotlib.font_manager import FontProperties

from Environment import Environment
from Export import Export


class Rail:
    def __init__(
        self,
        environment,
        rocket_mass,
        rocket_radius,
        aerodynamic_drag_coefficient,
        rail_length,
        rail_angle,
        thrust,
        frontal_area,
    ):
        self.environment = environment

        self.rocket_mass = rocket_mass
        self.rocket_radius = rocket_radius
        self.evaluate_frontal_area(frontal_area)
        self.aerodynamic_drag_coefficient = aerodynamic_drag_coefficient

        self.rail_length = rail_length
        self.rail_angle = rail_angle  # wrt. horizontal

        self.thrust = thrust

        # Initial differential equation conditions
        self.initial_position = 0.0
        self.initial_velocity = 0.0

        # Differential equation solution
        (self.time, self.position, self.velocity) = self.solve_rail()
        self.end_rail_velocity = self.solve_rail()[2][-1]

    def evaluate_frontal_area(self, frontal_area):
        """Frontal area computation from user input or cilindrical
        rocket approximation.

        Args:
            frontal_area (float): rocket cross sectional frontal area
        """
        if frontal_area is None:
            self.frontal_area = np.pi * (self.rocket_radius**2)
        else:
            self.frontal_area = frontal_area

    def vector_field(self, time, state_variables, parameters):
        """Generates the vector field of the corresponding simulations
        state variables (position from rail start, rocket velocity) as
        required for solve_ivp differential equation solver.

        Args:
            time (float): independent current time variable
            state_variables (list): simulation state variables
            to be solved
            parameters (tuple): main constants for vector state
            computation

        Returns:
            list: vector field of the state variables
        """
        position, velocity = state_variables
        (
            gravity,
            rocket_mass,
            frontal_area,
            aerodynamic_drag_coefficient,
            air_density,
            rail_angle,
            thrust,
        ) = parameters
        k_drag = (air_density * aerodynamic_drag_coefficient * frontal_area) / 2

        vector_state = [
            velocity,
            (thrust(time) - k_drag * velocity**2) / rocket_mass
            - gravity * np.sin(rail_angle),
        ]

        return vector_state

    def evaluate_jacobian(self, time, state_variables, parameters):
        """Explicit jacobian computation for faster computation
        of jacobian dependent solver methods, such as BDF or LSODA.

        Args:
            time (float): independent current time variable
            state_variables (list): simulation state variables
            to be solved
            parameters (tuple): main constants for vector state
            computation

        Returns:
            matrix: jacobian matrix
        """
        position, velocity = state_variables
        (
            gravity,
            rocket_mass,
            frontal_area,
            aerodynamic_drag_coefficient,
            air_density,
            rail_angle,
            thrust,
        ) = parameters
        k_drag = (air_density * aerodynamic_drag_coefficient * frontal_area) / 2

        dvelocity_dposition = 0
        dvelocity_dvelocity = 1
        dacceleration_dposition = 0
        dacceleration_dvelocity = -2 * k_drag * velocity / rocket_mass

        jacobian = [
            [dvelocity_dposition, dvelocity_dvelocity],
            [dacceleration_dposition, dacceleration_dvelocity],
        ]

        return jacobian

    def end_rail(self, time, state_variables, parameters):
        """Establishment of solver terminal conditions. The simulation ends
            if the rocket reaches the end of the rail, i.e. its position is
            greater than the rail length.

        Args:
            time (float): independent current time variable
            state_variables
            state_variables (list): simulation state variables
            to be solved
            parameters (tuple): main constants for vector state
        computation

        Returns:
            integer: boolean integer as termination parameter
        """
        position = state_variables[0]
        if position >= self.rail_length:
            return 0
        return 1

    end_rail.terminal = True

    def solve_rail(self):
        """Initial conditions setting and solver instatiation.

        Returns:
            object: solution object containing the solution
            for the differential equation state variables
        """

        state_variables = [self.initial_position, self.initial_velocity]

        parameters = [
            self.environment.gravity,
            self.rocket_mass,
            self.frontal_area,
            self.aerodynamic_drag_coefficient,
            self.environment.air_density,
            self.rail_angle,
            self.thrust,
        ]

        solution = solve_ivp(
            self.vector_field,
            (0.0, 100.0),
            state_variables,
            args=(parameters,),
            method="BDF",
            jac=self.evaluate_jacobian,
            events=self.end_rail,
            max_step=0.001,
        )

        self.solution = [solution.t, solution.y[0], solution.y[1]]

        return self.solution


class RailExport(Export):
    def __init__(self, Rail):
        self.Rail = Rail
        self.rail_exporting()

    def rail_exporting(self):
        """Method that calls Export class for solution exporting in a csv.

        Returns:
            None
        """
        try:
            Export.raw_simulation_data_export(
                self.Rail.solution,
                "data/rail_movement/rail_data.csv",
                ["Time", "Position", "Velocity"],
            )
        except OSError as err:
            print("OS error: {0}".format(err))
        return None

    def all_info(self):
        """Console logging of notable rail movement characteristics.

        Returns:
            None
        """
        print("Out rail velocity: {:.2f} m/s".format(self.Rail.end_rail_velocity))
        print("Mean rail velocity: {:.2f} m/s".format(np.mean(self.Rail.velocity)))
        print("Out rail time: {:.2f} s".format(self.Rail.time[-1]))

        return None

    def plotting(self):
        """Plot graphs of notable rail movement values.

        Returns:
            None
        """
        plt.figure(201, figsize=(16, 9))
        plt.plot(
            self.Rail.time,
            self.Rail.velocity,
            color="orange",
            linewidth=0.75,
            label=r"$v_{rail}$",
        )
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("velocity (m/s)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Rail velocity as function of time")
        plt.savefig("data/rail_movement/graphs/rail_velocity.png", dpi=200)

        plt.figure(202, figsize=(16, 9))
        plt.plot(
            self.Rail.time,
            self.Rail.position,
            color="orange",
            linewidth=0.75,
            label=r"$s_{rail}$",
        )
        plt.grid(True)
        plt.xlabel("time (s)")
        plt.ylabel("position (m)")
        plt.legend(prop=FontProperties(size=16))
        plt.title("Rail position as function of time")
        plt.savefig("data/rail_movement/graphs/rail_position.png", dpi=200)

        return None


if __name__ == "__main__":

    data_path = "data/burn_simulation/burn_data.csv"
    ext_time_list, ext_thrust = np.loadtxt(
        data_path, delimiter=",", unpack=True, skiprows=1, usecols=(0, 4)
    )

    thrust = interpolate.interp1d(ext_time_list, ext_thrust)

    Pirassununga = Environment(-0.38390456, 627, ellipsoidal_model=True)
    Keron_test = Rail(Pirassununga, 10, 100 / 2000, 0.5, 4, np.pi / 2, thrust, None)

    RailExport(Keron_test).all_info()
    RailExport(Keron_test).plotting()
