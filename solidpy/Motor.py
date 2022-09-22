# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = "MIT"
_license_ = ""

import numpy as np


class Motor:
    def __init__(
        self,
        grain,
        grain_number,
        chamber_inner_radius,
        nozzle_throat_radius,
        nozzle_exit_radius,
        nozzle_angle,
        chamber_length=None,
    ):
        self.grain = grain
        self.grain_number = grain_number
        self.evaluate_chamber_length(chamber_length)
        self.chamber_area = np.pi * chamber_inner_radius**2
        self.nozzle_throat_area = np.pi * nozzle_throat_radius**2
        self.nozzle_exit_area = np.pi * nozzle_exit_radius**2
        self.nozzle_angle = nozzle_angle
        self.expansion_ratio = self.nozzle_exit_area / self.nozzle_throat_area
        self.evaluate_chamber_volume()
        self.evaluate_propellant_volume()
        self.evaluate_free_volume()
        self.evaluate_total_burn_area()
        self.evaluate_Kn()

    def evaluate_chamber_length(self, chamber_length):
        if chamber_length is None:
            self.chamber_length = self.grain.initial_height * self.grain_number
        else:
            self.chamber_length = chamber_length

    def evaluate_chamber_volume(self):
        self.chamber_volume = self.chamber_area * self.chamber_length
        return self.chamber_volume

    def evaluate_propellant_volume(self):
        self.propellant_volume = self.grain_number * self.grain.volume
        return self.propellant_volume

    def evaluate_free_volume(self):
        self.free_volume = (
            self.evaluate_chamber_volume() - self.evaluate_propellant_volume()
        )
        return self.free_volume

    def evaluate_Kn(self):
        self.Kn = (self.grain_number * self.grain.burn_area) / self.nozzle_throat_area
        return self.Kn

    def evaluate_total_burn_area(self):
        self.total_burn_area = self.grain_number * self.grain.burn_area
