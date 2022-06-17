# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import csv
import numpy as np
from scipy import interpolate
from scipy.integrate import solve_ivp

from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, grid, legend, title, savefig, show

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
        if frontal_area is None:
            self.frontal_area = np.pi * (self.rocket_radius**2)
        else:
            self.frontal_area = frontal_area

    def vector_field(self, time, state_variables, parameters):
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

    def solve_rail(self):

        """Variables"""
        state_variables = [self.initial_position, self.initial_velocity]

        """Parameters"""
        parameters = [
            self.environment.gravity,
            self.rocket_mass,
            self.frontal_area,
            self.aerodynamic_drag_coefficient,
            self.environment.air_density,
            self.rail_angle,
            self.thrust,
        ]

        def end_rail(time, state_variables, parameters):
            position = state_variables[0]
            if position >= self.rail_length:
                return 0
            return 1

        end_rail.terminal = True

        """Solver"""
        solution = solve_ivp(
            self.vector_field,
            (0.0, 100.0),
            state_variables,
            args=(parameters,),
            method="BDF",
            jac=self.evaluate_jacobian,
            events=end_rail,
            max_step=0.001,
        )

        self.solution = [solution.t, solution.y[0], solution.y[1]]

        return self.solution


class RailExport(Export):
    def __init__(self, Rail):
        self.Rail = Rail
        self.post_processing()

    def post_processing(self):
        Export.raw_simulation_data_export(
            self.Rail.solution,
            "data/rail_movement/rail_data.csv",
            ["Time", "Position", "Velocity"],
        )
        return None

    def all_info(self):
        print("Out rail velocity: {:.2f} m/s".format(self.Rail.end_rail_velocity))
        print("Mean rail velocity: {:.2f} m/s".format(np.mean(self.Rail.velocity)))
        print("Out rail time: {:.2f} s".format(self.Rail.time[-1]))

        return None

    def plotting(self):

        figure(1, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            self.Rail.time, self.Rail.velocity, "b", linewidth=0.75, label=r"$v_{rail}$"
        )
        legend(prop=FontProperties(size=16))
        title("Rail velocity as function of time")
        savefig("data/rail_movement/graphs/rail_velocity.png", dpi=200)

        figure(2, figsize=(16, 9))
        xlabel("t")
        grid(True)
        plot(
            self.Rail.time, self.Rail.position, "b", linewidth=0.75, label=r"$s_{rail}$"
        )
        legend(prop=FontProperties(size=16))
        title("Rail position as function of time")
        savefig("data/rail_movement/graphs/rail_position.png", dpi=200)

        return None


if __name__ == "__main__":

    data_path = "data/burn_simulation/burn_data.csv"
    ext_time_list, ext_thrust = np.loadtxt(
        data_path, delimiter=",", unpack=True, skiprows=1, usecols=(0, 4)
    )

    thrust = interpolate.interp1d(ext_time_list, ext_thrust)

    Pirassununga = Environment(-0.38390456, 750, ellipsoidal_model=True)
    Keron_test = Rail(Pirassununga, 10, 100 / 2000, 0.5, 4, np.pi / 2, thrust, None)

    RailExport(Keron_test).all_info()
    RailExport(Keron_test).plotting()
