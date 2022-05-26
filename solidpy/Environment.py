# -*- coding: utf-8 -*-

import math
import scipy.constants as const


class Environment:
    def __init__(
        self,
        atmospheric_pressure=const.atmosphere,
        air_density=1.225,
        latitude=0,
        gravity=None,
        ellipsoidal_model=False,
    ):
        self.atmospheric_pressure = atmospheric_pressure
        self.air_density = air_density
        self.latitude = latitude
        self.gravity = self.evaluate_gravity(gravity, latitude, ellipsoidal_model)

    def evaluate_gravity(self, gravity, latitude, ellipsoidal_model):
        if gravity is None and ellipsoidal_model:
            # Somigliana approximation for ellipsoidal gravity
            gravity = 9.78032533590389 * (
                (1 + 0.001931852652458 * (math.sin(latitude)) ** 2)
                / (math.sqrt(1 - 0.006694379990141 * (math.sin(latitude)) ** 2))
            )
            return gravity
        elif not ellipsoidal_model:
            return const.g
        return gravity
