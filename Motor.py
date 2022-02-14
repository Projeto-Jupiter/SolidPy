# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import numpy as np
from Grain import Grain


class Motor:
    def __init__(
        self,
        grain,
        grain_number,
        chamber_inner_radius,
        nozzle_throat_radius,
        chamber_length=None,
        nozzle_exit_radius=None,
    ):
        self.grain = grain
        self.grain_number = grain_number
        self.evaluate_chamber_length(chamber_length)
        self.chamber_area = np.pi * chamber_inner_radius**2
        self.nozzle_throat_area = np.pi * nozzle_throat_radius**2
        # self.nozzle_exit_radius = np.pi*nozzle_exit_radius**2
        self.evaluate_chamber_volume()
        self.evaluate_propellant_volume()
        self.evaluate_free_volume()
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


""" Grao_Leviata = Grain(outer_radius=71.92 / 2000, initial_inner_radius=31.92 / 2000)
Leviata = Motor(
    Grao_Leviata,
    grain_number=4,
    chamber_inner_radius=77.92 / 2000,
    nozzle_throat_radius=8.75 / 2000,
)
print(Leviata.chamber_volume)
print(Leviata.propellant_volume)
print(Leviata.free_volume)
print(Leviata.Kn)
 """
