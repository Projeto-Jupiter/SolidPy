# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, João Lemes Gribel Soares, Thais Silva Melo, Tiago Mariotto Lucio"
_copyright_ = "x"
_license_ = "x"

import numpy as np
from abc import ABC, abstractmethod


class Grain(ABC):
    def __init__(
        self, outer_radius, inner_radius, height=None, regressed_length=0, mass=None
    ):
        self.outer_radius = outer_radius
        self.inner_radius = inner_radius
        self.initial_inner_radius = self.inner_radius
        self.height = height
        self.initial_height = self.height
        self.initial_volume = self.volume
        self.mass = mass
        self.regressed_length = regressed_length

    @property
    def inner_radius(self):
        return self._inner_radius

    @inner_radius.setter
    def inner_radius(self, radius):
        if radius < self.outer_radius:
            self._inner_radius = radius
        else:
            raise Exception("Inner radius cannot be greater than outer radius.")

    @property
    def height(self):
        return self._height

    @height.setter
    @abstractmethod
    def height(self, height):
        pass

    @property
    def regressed_length(self):
        return self._regressed_length

    @regressed_length.setter
    def regressed_length(self, length):
        self._regressed_length = length
        self.burn_regress(length)

    @abstractmethod
    def burn_regress(self, regression):
        pass

    @property
    @abstractmethod
    def contour_length(self):
        pass

    @property
    @abstractmethod
    def transversal_area(self):
        pass

    @property
    def burn_area(self):
        return 2 * self.transversal_area + self.contour_length * self.height

    @property
    def volume(self):
        return self.transversal_area * self.height


class Bates(Grain):
    def __init__(
        self, outer_radius, inner_radius, height=None, regressed_length=0, mass=None
    ):
        super().__init__(outer_radius, inner_radius, height, regressed_length, mass)

    @Grain.height.setter
    def height(self, height):
        if height is None:
            self._height = 3 * self.outer_radius + self.inner_radius
        else:
            self._height = height

    def burn_regress(self, length):
        self._inner_radius = self.initial_inner_radius + length
        self._height = self.initial_height - 2 * length

    @Grain.contour_length.getter
    def contour_length(self):
        return 2 * np.pi * self.inner_radius

    @Grain.transversal_area.getter
    def transversal_area(self):
        return np.pi * (self.outer_radius**2 - self.inner_radius**2)


class CustomGeometry(Grain):
    def __init__(self, contour_mesh):
        ...

    def set_parametric_contour(self):
        ...

    def evaluate_contour_length(self):
        ...

    def evaluate_transversal_area(self):
        ...


if __name__ == "__main__":
    Grao_Leviata = Bates(outer_radius=71.92 / 2000, inner_radius=31.92 / 2000)
    print(Grao_Leviata.burn_area)
    Grao_Leviata.burn_regress((71.92 / 2000 - 31.92 / 2000) / 2)
    print(Grao_Leviata.burn_area)
