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
        self.transversal_area = np.pi * chamber_inner_radius**2
        self.chamber_length = chamber_length
        self.nozzle_throat_area = np.pi * nozzle_throat_radius**2
        self.nozzle_exit_area = np.pi * nozzle_exit_radius**2
        self.nozzle_angle = nozzle_angle
        self.expansion_ratio = self.nozzle_exit_area / self.nozzle_throat_area

    @property
    def chamber_length(self):
        return self._chamber_length

    @chamber_length.setter
    def chamber_length(self, length):
        if length is None:
            self._chamber_length = self.grain.height * self.grain_number
        else:
            self._chamber_length = length

    @property
    def chamber_volume(self):
        return self.transversal_area * self.chamber_length

    @property
    def propellant_volume(self):
        return self.grain_number * self.grain.volume

    @property
    def free_volume(self):
        return self.chamber_volume - self.propellant_volume

    @property
    def Kn(self):
        self.Kn = (self.grain_number * self.grain.burn_area) / self.nozzle_throat_area
        return self.Kn

    @property
    def total_burn_area(self):
        return self.grain_number * self.grain.burn_area
