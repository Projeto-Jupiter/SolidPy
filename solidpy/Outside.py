# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, Luís Felipe Biancardi Palharini, Matheus Larrondo Portiolli, Pedro Henrique Marinho Bressan"
_copyright_ = "MIT"
_license_ = "x"

import csv
import numpy as np
from scipy import interpolate
from scipy.integrate import solve_ivp


class Outside:
    def __init__(
        self,
        atmospheric_pressure,
        rocket_mass,
        rocket_radius,
        aerodynamic_drag_coefficient,
        air_density,
        rail_length,
        rail_angle,
        frontal_area,
        thrust,
        gravity=9.8,
    ):

        # Constant parameters
        self.atmospheric_pressure = atmospheric_pressure
        self.gravity = gravity
        self.rocket_mass = rocket_mass
        self.rocket_radius = rocket_radius
        self.frontal_area = None
        self.aerodynamic_drag_coefficient = aerodynamic_drag_coefficient
        self.air_density = air_density
        self.rail_length = rail_length
        self.rail_angle = rail_angle  # wrt. horizontal
        self.evaluate_frontal_area(frontal_area)

        # Thrust vector temporary state
        # pending future "Burn" class
        self.thrust = thrust

        # Initial differential equation conditions
        self.initial_position = 0.0
        self.initial_velocity = 0.0
        self.time_span = 0.0, 2.5
        # M for Manual, A for Auto, ESn for Evenly Spaced (n=number steps)
        self.time_method = "M"

        # Differential equation solution
        self.end_rail_velocity = self.solve_rail()

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
            (thrust(time) - k_drag * velocity**2) / (rocket_mass)
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

    def set_time_steps(self):
        time_method = self.time_method
        time_span = self.time_span
        if time_method == "M":
            return ext_time_list
        elif time_method[:2] == "ES":
            number_steps = int(time_method[2:])
            time_steps = np.linspace(time_span[0], time_span[1], number_steps)
            return time_steps
        return None

    def solve_rail(self):

        """Variables"""
        state_variables = [self.initial_position, self.initial_velocity]

        """Parameters"""
        parameters = [
            self.gravity,
            self.rocket_mass,
            self.frontal_area,
            self.aerodynamic_drag_coefficient,
            self.air_density,
            self.rail_angle,
            self.thrust,
        ]
        absolute_error = 1.0e-12
        relative_error = 1.0e-8

        def end_rail(time, state_variables, parameters):
            position = state_variables[0]
            if position >= self.rail_length:
                return 0
            return 1

        end_rail.terminal = True

        """Solver"""
        # uncomment "atol" and "rtol" for manual error control
        solution = solve_ivp(
            self.vector_field,
            self.time_span,
            state_variables,
            method="BDF",
            t_eval=self.set_time_steps(),
            args=(parameters,),
            events=end_rail,
            max_step=0.005,
            jac=self.evaluate_jacobian
            # atol=absolute_error,
            # rtol=relative_error,
        )

        """Export data"""
        with open("data/rail_movement/rail_data.csv", "w") as rail_data:
            solution_writer = csv.writer(rail_data)
            solution_writer.writerow(["Time", "Position", "Velocity"])
            for solution_time, solution_position, solution_velocity in zip(
                solution.t, solution.y[0], solution.y[1]
            ):
                solution_writer.writerow(
                    [solution_time, solution_position, solution_velocity]
                )

        return solution_velocity


"""Simulation for testing and debugging purposes only"""

# data_path = "data/static_fires/Keron.csv"
# ext_time_list, ext_thrust = np.loadtxt(
#    data_path, delimiter=";", unpack=True, skiprows=1
# )

data_path = "data/burn_simulation/burn_data.csv"
ext_time_list, ext_thrust = np.loadtxt(
    data_path, delimiter=",", unpack=True, skiprows=1, usecols=(0, 1)
)

thrust = interpolate.interp1d(ext_time_list, ext_thrust)

Keron_test = Outside(10e5, 10, 21.4e-3, 0.5, 1.25, 4, np.pi / 2, None, thrust)
print(Keron_test.end_rail_velocity)
