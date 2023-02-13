# -*- coding: utf-8 -*-

_author_ = "Caio Eduardo dos Santos de Souza, João Lemes Gribel Soares, Pedro Bressan, Thais Silva Melo, Tiago Mariotto Lucio"
_copyright_ = "x"
_license_ = "x"

import skfmm
import numpy as np
import pylab as plt
import scipy.interpolate as interpolate

from abc import ABC, abstractmethod
from matplotlib import animation

from .Export import Export


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
        if length:
            self.burn_regress(length)

    def burn_regress(self, regressed_length):
        self._inner_radius = self.initial_inner_radius + regressed_length
        self._height = self.initial_height - 2 * regressed_length

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

    @Grain.contour_length.getter
    def contour_length(self):
        return 2 * np.pi * self.inner_radius

    @Grain.transversal_area.getter
    def transversal_area(self):
        return np.pi * (self.outer_radius**2 - self.inner_radius**2)


class MarchingGrain(Grain):
    def __init__(
        self,
        outer_radius,
        height=None,
        regressed_length=0,
        mass=None,
        contour_number=25,
        grid_refinement=1000,
    ):
        super().__init__(outer_radius, height, regressed_length, mass)
        self.contour_number = contour_number
        self.grid_refinement = grid_refinement

        self.x_mesh, self.y_mesh = self.generate_grid()
        self.mesh_radius = np.sqrt(self.x_mesh**2 + self.y_mesh**2)

        self.inner_contour = self.evaluate_regression()
        (
            self.total_contour_length,
            self.total_contour_area,
        ) = self.evaluate_contour_properties()

    def generate_grid(self):
        partition, self.grid_step = np.linspace(
            -self.outer_radius,
            self.outer_radius,
            self.grid_refinement + 1,
            retstep=True,
        )
        return np.meshgrid(partition, partition)

    def bound_regression(self, contour):
        limit_outer_contour = self.mesh_radius >= self.outer_radius
        limit_inner_contour = 2 * self.grid_step < -contour

        contour = np.ma.MaskedArray(contour, limit_outer_contour)
        contour = np.ma.MaskedArray(contour, limit_inner_contour)
        return contour

    def evaluate_regression(self):
        # Create inner contour and bound its regression by masking
        inner_contour = self.generate_inner_shape()
        inner_contour = self.bound_regression(inner_contour)

        # Calculate regression from contour
        distance_grid = skfmm.distance(inner_contour, float(self.grid_step))
        regression_contours = plt.contour(
            self.x_mesh, self.y_mesh, distance_grid, self.contour_number
        )
        return regression_contours

    def evaluate_contour_properties(self):
        length = []
        area = []
        regression_steps = np.linspace(0, self.outer_radius, self.contour_number)

        for contour in self.inner_contour.collections[1:-1]:
            x = contour.get_paths()[0].vertices[:, 0]
            y = contour.get_paths()[0].vertices[:, 1]
            # Sum line elements legths by Pythagoras
            length.append(np.sum((np.diff(x) ** 2 + np.diff(y) ** 2) ** (1 / 2)))
            # Evaluate area from points by Green's Theorem
            area.append(np.abs(0.5 * np.sum(y[:-1] * np.diff(x) - x[:-1] * np.diff(y))))

        self.regression_steps = regression_steps[: min(len(length), len(area))]

        return (
            interpolate.interp1d(self.regression_steps, length, kind="cubic"),
            interpolate.interp1d(self.regression_steps, area, kind="cubic"),
        )

    @Grain.contour_length.getter
    def contour_length(self):
        return self.total_contour_length(self.regressed_length)

    @Grain.transversal_area.getter
    def transversal_area(self):
        return np.pi * self.outer_radius**2 - self.total_contour_area(
            self.regressed_length
        )

    @abstractmethod
    def generate_inner_shape(self):
        pass


class Star(MarchingGrain):
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
        grid_refinement=1000,
    ):
        self.star_maximum = star_maximum
        self.star_minimum = star_minimum
        self.star_points = star_points
        self.inner_radius = star_minimum
        self.initial_inner_radius = star_minimum

        super().__init__(
            outer_radius,
            height,
            regressed_length,
            mass,
            contour_number,
            grid_refinement,
        )

    @Grain.height.setter
    def height(self, height):
        if height is None:
            self._height = (
                3 * self.outer_radius + (self.star_maximum + self.star_minimum) / 2
            )
        else:
            self._height = height

    def generate_inner_shape(self):
        angle = np.arctan2(self.y_mesh, self.x_mesh)
        standard_star_contour = (
            self.mesh_radius
            - (self.star_maximum + self.star_minimum) / 2
            - (self.star_maximum - self.star_minimum)
            / 2
            * np.sin(self.star_points * angle)
        )
        return standard_star_contour


class CustomGeometry(MarchingGrain):
    def __init__(
        self,
        outer_radius,
        inner_radius,
        height,
        input_method="polar",
        regressed_length=0,
        mass=None,
        contour_number=25,
        grid_refinement=1000,
    ):
        self.input_method = input_method
        self.inner_radius = inner_radius

        self.method_map = {
            "polar": self.generate_polar,
        }

        super().__init__(
            outer_radius,
            height,
            regressed_length,
            mass,
            contour_number,
            grid_refinement,
        )

    @Grain.height.setter
    def height(self, height):
        self._height = height

    def generate_polar(self):
        mesh_angle = np.arctan2(self.y_mesh, self.x_mesh)
        input_radius = np.abs(self.inner_radius(mesh_angle))
        self.initial_inner_radius = np.min(input_radius)
        return self.mesh_radius - input_radius

    def generate_parametric(self):
        # radius = np.sqrt(self.x_mesh**2 + self.y_mesh**2)
        # x, y = self.inner_radius(np.linspace())
        Ellipsis

    def generate_mesh(self):
        x, y = self.inner_radius[:, 0], self.inner_radius[:, 1]
        input_radius = np.sqrt(x**2 + y**2)
        self.initial_inner_radius = np.min(np.sqrt(x**2 + y**2))
        return self.mesh_radius - input_radius

    def generate_inner_shape(self):
        try:
            return self.method_map.get(self.input_method)()
        except ValueError:
            raise ("Check inner_radius input.")


class GrainExport(Export):
    def __init__(self, grain):
        self.grain = grain

    def plotting(self):
        plt.figure(301, figsize=(16, 9))
        plt.contour(self.grain.inner_contour)
        plt.gca().set_aspect(1)
        plt.title("Grain regression from its zero level port")
        plt.show()
        plt.close()

        plt.figure(302, figsize=(16, 9))
        plt.plot(
            self.grain.regression_steps,
            self.grain.total_contour_length(self.grain.regression_steps),
        )
        plt.title("Contour length as a function of regression")
        plt.show()
        plt.close()

        plt.figure(303, figsize=(16, 9))
        plt.plot(
            self.grain.regression_steps,
            self.grain.total_contour_area(self.grain.regression_steps),
        )
        plt.title("Transversal area as a function of regression")
        plt.show()
        plt.close()

    def animation(self, index):
        contour = self.grain.inner_contour.collections[1:-1][index]
        x = contour.get_paths()[0].vertices[:, 0]
        y = contour.get_paths()[0].vertices[:, 1]
        return plt.plot(x, y)

    def export_animation(self):
        figure = plt.figure()
        anim = animation.FuncAnimation(
            figure,
            self.animation,
            frames=len(self.grain.regression_steps),
            interval=2000,
        )
        anim.show()


if __name__ == "__main__":
    Grao_Leviata = Bates(outer_radius=71.92 / 2000, inner_radius=31.92 / 2000)
    print(Grao_Leviata.burn_area)
    Grao_Leviata.burn_regress((71.92 / 2000 - 31.92 / 2000) / 2)
    print(Grao_Leviata.burn_area)

    Star_Test = Star(
        outer_radius=30 / 1000, star_maximum=50 / 3000, star_minimum=10 / 3000
    )
