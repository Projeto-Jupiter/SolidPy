# -*- coding: utf-8 -*-

import math
import scipy.constants as const
from ambiance import Atmosphere


class Environment:
    def __init__(
        self,
        latitude=0,
        altitude=0,
        gravity=None,
        ellipsoidal_model=False,
    ):
        atmosphere = Atmosphere(altitude)
        self.atmospheric_pressure = atmosphere.pressure[0]
        self.air_density = atmosphere.density[0]
        self.latitude = latitude
        self.standard_gravity = const.g
        self.gravity = self.evaluate_gravity(
            gravity, latitude, altitude, ellipsoidal_model
        )

    def evaluate_gravity(self, gravity, latitude, altitude, ellipsoidal_model):
        """Local gravity acceleration evaluation based on given user input.

        Args:
            gravity (float): optional user input gravity value
            latitude (float): local latitude
            altitude (float): local altitude
            ellipsoidal_model (bool): user input whether an ellipsoidal
            Earth ought be simulated

        Returns:
            float: local gravity acceleration
        """
        if gravity is None and ellipsoidal_model:
            # Somigliana approximation for ellipsoidal gravity
            # Source: "National Geospatial-intelligence Agency
            # Standardization Document" "NGA.STND.0036_1.0.0_WGS84"
            a = 6378137.0  # semi_major_axis
            f = 1 / 298.257223563  # flattening_factor
            m_rot = 3.449786506841e-3  # rotation_factor
            g_e = 9.7803253359  # normal gravity at equator
            k_somgl = 1.931852652458e-3  # normal gravity formula const.
            first_ecc_sqrd = 6.694379990141e-3  # square of first eccentricity

            gravity_somgl = g_e * (
                (1 + k_somgl * (math.sin(latitude)) ** 2)
                / (math.sqrt(1 - first_ecc_sqrd * (math.sin(latitude)) ** 2))
            )
            height_correction = (
                1
                - 2 / a * (1 + f + m_rot - 2 * f * (math.sin(latitude)) ** 2) * altitude
                + 3 * altitude**2 / a**2
            )

            return gravity_somgl * height_correction
        elif not ellipsoidal_model:
            return self.standard_gravity
        return gravity
