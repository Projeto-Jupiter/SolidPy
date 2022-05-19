# -*- coding: utf-8 -*-

import math


class Environment:
    def __init__(self, atmospheric_pressure, air_density, latitude=0, gravity=None):
        self.atmospheric_pressure = atmospheric_pressure
        self.air_density = air_density
        self.latitude = latitude
        self.gravity = self.evaluate_gravity(gravity)

    def evaluate_gravity(self, gravity):
        if gravity is None:
            # Somigliana approximation for ellipsoidal gravity
            gravity = 9.78032533590389 * (
                (1 + 0.001931852652458 * (math.sin(self.latitude)) ** 2)
                / (math.sqrt(1 - 0.006694379990141 * (math.sin(self.latitude)) ** 2))
            )
            return gravity
        return gravity
