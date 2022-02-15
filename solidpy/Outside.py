# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, Luís Felipe Biancardi Palharini, Matheus Larrondo Portiolli, Pedro Henrique Marinho Bressan"
_copyright_ = "x"
_license_ = "x"


import numpy as np
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
        self.evaluate_frontal_area(self, frontal_area)

        # Thrust vector temporary state
        # pending future "Burn" class
        self.thrust = thrust

        # Initial state variable conditions
        self.initial_position = 0.0
        self.initial_velocity = 0.0

        # Differential equation solution
        self.end_rail_velocity = self.solve_rail()

    def evaluate_frontal_area(self, frontal_area):
        if frontal_area is None:
            frontal_area = np.pi * (self.rocket_radius**2)
        else:
            self.frontal_area = frontal_area

    def vector_field(self, state_variables, time, parameters):
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
            (thrust - k_drag * velocity**2) / (rocket_mass)
            - gravity * np.sin(rail_angle),
        ]

        return vector_state

    def solve_rail(self):
        state_variables = [self.initial_position, self.initial_velocity]
        parameters = [
            self.gravity,
            self.rocket_mass,
            self.frontal_area,
            self.aerodynamic_drag_coefficient,
            self.air_density,
            self.rail_angle,
            self.thrust,
        ]
        time_span = 0.0, 1.0  # initial,  final
        absolute_error = 1.0e-12
        relative_error = 1.0e-8

        solution = solve_ivp(
            self.vector_field,
            time_span,
            state_variables,
            method="RK45",
            args=(parameters,),
            # atol=absolute_error, uncomment these lines
            # rtol=relative_error, for manual error control
        )

        with open("rail_data.dat", "w") as rail_data:
            for solution_var, solution_time in zip(solution.y, solution.t):
                if solution_var[0] > self.rail_length:
                    break
                else:
                    print(
                        solution_time, solution_var[0], solution_var[1], file=rail_data
                    )

        return solution_var[1]
