# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, João Lemes Gribel Soares, Thais Silva Melo"
_copyright_ = "x"
_license_ = "x"

import numpy as np


class Grain:
    def __init__(
        self,
        outer_radius,
        initial_inner_radius,
        geometry="tubular",
        initial_height=None,
    ):

        self.outer_radius = outer_radius
        self.initial_inner_radius = initial_inner_radius
        self.evaluate_grain_height(initial_height)
        self.geometry = geometry
        self.evaluate_grain_geometry()

    def evaluate_grain_height(self, initial_height):
        if initial_height is None:
            self.initial_height = 3 * self.outer_radius + self.initial_inner_radius  ### Verificar
        else:
            self.initial_height = initial_height

        return self.initial_height

    def evaluate_grain_geometry(self):
        if self.geometry == "tubular":
            self.evaluate_tubular_grain_area()
        elif self.geometry == "star":
            pass
        else:
            print("Not a valid geometry type")

    def evaluate_tubular_grain_area(self):
        transversal_area = 2 * np.pi * (self.outer_radius ** 2 - self.initial_inner_radius ** 2)
        longitudinal_area = 2 * np.pi * self.initial_inner_radius * self.initial_height
        self.grain_area = transversal_area + longitudinal_area
        return self.grain_area
