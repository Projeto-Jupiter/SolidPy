# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, João Lemes Gribel Soares, Thais Silva Melo, Tiago Mariotto Lucio"
_copyright_ = "x"
_license_ = "x"

import numpy as np


class Grain:
    def __init__(
        self,
        outer_radius,
        initial_inner_radius,
        initial_height=None,
        geometry="tubular",
    ):

        self.outer_radius = outer_radius
        self.inner_radius = initial_inner_radius
        self.evaluate_grain_initial_height(initial_height)
        self.height = self.initial_height
        self.geometry = geometry
        self.evaluate_grain_geometry()
        self.evaluate_grain_volume()

    def evaluate_grain_initial_height(self, initial_height):
        if initial_height is None:
            self.initial_height = 3 * self.outer_radius + self.inner_radius
        else:
            self.initial_height = initial_height

    def evaluate_grain_geometry(self):
        if self.geometry == "tubular":
            self.evaluate_tubular_burn_area()
        elif self.geometry == "star":
            pass
        else:
            print("Not a valid geometry type")

    def evaluate_tubular_burn_area(self):
        transversal_area = 2 * np.pi * (self.outer_radius**2 - self.inner_radius**2)
        longitudinal_area = 2 * np.pi * self.inner_radius * self.height
        self.burn_area = transversal_area + longitudinal_area

        return self.burn_area

    def evaluate_grain_volume(self):
        self.volume = (
            np.pi * (self.outer_radius**2 - self.inner_radius**2) * self.height
        )
        return self.volume
