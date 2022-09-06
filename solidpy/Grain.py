# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, João Lemes Gribel Soares, Pedro Bressan, Thais Silva Melo, Tiago Mariotto Lucio"
_copyright_ = "x"
_license_ = "x"

import skfmm
import numpy as np
import pylab as plt

# from typing import Callable
from abc import ABC, abstractmethod
import scipy.interpolate as interpolate


class Grain(ABC):
    def __init__(self, outer_radius, height=None, regressed_length=0, mass=None):
        self.outer_radius = outer_radius
        self.height = height
        self.initial_height = self.height
        self.mass = mass
        self.regressed_length = regressed_length

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
    def burn_regress(self, regressed_length, burnrate):
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
        self,
        outer_radius,
        inner_radius,
        height=None,
        regressed_length=0,
        mass=None,
    ):
        self.inner_radius = inner_radius
        self.initial_inner_radius = self.inner_radius
        super().__init__(outer_radius, height, regressed_length, mass)

    @property
    def inner_radius(self):
        return self._inner_radius

    @inner_radius.setter
    def inner_radius(self, radius):
        self._inner_radius = radius

    @Grain.height.setter
    def height(self, height):
        if height is None:
            self._height = 3 * self.outer_radius + self.inner_radius
        else:
            self._height = height

    def burn_regress(self, regressed_length):
        self._inner_radius = self.initial_inner_radius + regressed_length
        self._height = self.initial_height - 2 * regressed_length

    @Grain.contour_length.getter
    def contour_length(self):
        return 2 * np.pi * self.inner_radius

    @Grain.transversal_area.getter
    def transversal_area(self):
        return np.pi * (self.outer_radius**2 - self.inner_radius**2)


class Star(Grain):
    def __init__(
        self,
        outer_radius,
        star_maximum,
        star_minimum,
        star_points=5,
        height=None,
        regressed_length=0,
        mass=None,
        contour_number=25,
        grid_refinement=1001,
    ):

        self.star_maximum = star_maximum
        self.star_minimum = star_minimum
        self.star_points = star_points
        self.contour_number = contour_number
        self.grid_refinement = grid_refinement
        self.inner_radius = star_minimum
        self.initial_inner_radius = star_minimum
        super().__init__(outer_radius, height, regressed_length, mass)

        self.inner_contour = self.set_inner_contour()
        (
            self.total_contour_length,
            self.total_contour_area,
        ) = self.evaluate_contour_properties()

    @Grain.height.setter
    def height(self, height):
        if height is None:
            self._height = (
                3 * self.outer_radius + (self.star_maximum + self.star_minimum) / 2
            )
        else:
            self._height = height

    def set_inner_contour(self):

        grid = np.linspace(
            -self.outer_radius, self.outer_radius, self.grid_refinement + 1
        ), np.linspace(-self.outer_radius, self.outer_radius, self.grid_refinement + 1)
        x_mesh, y_mesh = np.meshgrid(*grid)

        polar_radius = (x_mesh**2 + y_mesh**2) ** (1 / 2)
        polar_angle = np.arctan2(y_mesh, x_mesh)

        minimum_step = 2 * self.outer_radius / self.grid_refinement

        inner_radius_zero_contour = (
            polar_radius
            - (self.star_maximum + self.star_minimum) / 2
            - (self.star_maximum - self.star_minimum)
            / 2
            * np.sin(self.star_points * polar_angle)
        )

        limit_outer_contour = x_mesh**2 + y_mesh**2 >= self.outer_radius**2
        limit_inner_contour = (
            polar_radius
            < (self.star_maximum + self.star_minimum) / 2
            + (self.star_maximum - self.star_minimum)
            / 2
            * np.sin(self.star_points * polar_angle)
            - 10 * minimum_step
        )

        inner_radius_zero_contour = np.ma.MaskedArray(
            inner_radius_zero_contour, limit_outer_contour
        )
        inner_radius_zero_contour = np.ma.MaskedArray(
            inner_radius_zero_contour, limit_inner_contour
        )

        distance_grid = skfmm.distance(inner_radius_zero_contour, minimum_step)

        contours = plt.contour(x_mesh, y_mesh, distance_grid, self.contour_number)

        return contours

    def burn_regress(self, regressed_length):
        self.inner_radius = self.initial_inner_radius + regressed_length
        self._height = self.initial_height - 2 * regressed_length

    def evaluate_contour_properties(self):
        length = []
        area = []
        regression_steps = np.linspace(0, self.outer_radius, self.contour_number)

        for contour in self.inner_contour.collections[1:-1]:
            x = contour.get_paths()[0].vertices[:, 0]
            y = contour.get_paths()[0].vertices[:, 1]
            length.append(np.sum((np.diff(x) ** 2 + np.diff(y) ** 2) ** (1 / 2)))
            area.append(np.abs(0.5 * np.sum(y[:-1] * np.diff(x) - x[:-1] * np.diff(y))))

        return (
            interpolate.interp1d(regression_steps[: len(length)], length, kind="cubic"),
            interpolate.interp1d(regression_steps[: len(area)], area, kind="cubic"),
        )

    @Grain.contour_length.getter
    def contour_length(self):
        return self.total_contour_length(self.regressed_length)

    @Grain.transversal_area.getter
    def transversal_area(self):
        return np.pi * self.outer_radius**2 - self.total_contour_area(
            self.regressed_length
        )


class CustomGeometry(Grain):
    def __init__(
        self, outer_radius, inner_radius, height=None, regressed_length=0, mass=None
    ):
        super().__init__(outer_radius, inner_radius, height, regressed_length, mass)
        ...

    def set_inner_contour(self):
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

    Star_Test = Star(
        outer_radius=30 / 1000, star_maximum=50 / 3000, star_minimum=10 / 3000
    )
