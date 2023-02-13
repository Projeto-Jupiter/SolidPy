# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import math

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties

from .Grain import Bates, Star
from .Motor import Motor
from .Environment import Environment
from .Propellant import Propellant
from .Burn import Burn
from .Export import Export


class BurnEmpirical(Burn):
    def __init__(
        self,
        motor,
        propellant,
        environment=Environment(),
        empirical_data=None,
    ):
        Burn.__init__(self, motor, propellant, environment)
        self.grain = motor.grain
        self.empirical_propellant_density = self.evaluate_propellant_density()

        """Empirical known results (e.g. static-fire)"""
        self.empirical_time_steps, self.empirical_thrust = empirical_data
        self.empirical_chamber_pressure = self.evaluate_empirical_chamber_pressure()
        self.empirical_burn_rate = self.evaluate_empirical_burn_rate()

    def evaluate_propellant_density(self):
        """Method that overwrites propellant standard density for user
        empirical evaluation.

        Returns:
            float: propellant density
        """
        if self.grain.mass is not None:
            return self.grain.mass / self.grain.volume
        return self.propellant.density

    def evaluate_empirical_chamber_pressure(self):
        """Computation of motor's chamber pressure from user supplied
        thrust data.

        Source:
        https://www.grc.nasa.gov/www/k-12/rocket/rktthsum.html

        Returns:
            list: chamber pressure values for each thrust input
        """
        T_0, R, _, k, A_t = self.parameters
        chamber_pressure_list = []

        for thrust in self.empirical_thrust:
            current_chamber_pressure = (
                thrust + self.motor.nozzle_exit_area * self.environment_pressure
            ) / (
                A_t
                * np.sqrt(k / (R * T_0))
                * math.pow((2 / (k + 1)), ((k + 1) / (2 * (k - 1))))
                * self.evaluate_exit_velocity()
                + math.pow(
                    (1 + (k - 1) / 2 * self.evaluate_exit_mach() ** 2), -k / (k - 1)
                )
                * self.motor.nozzle_exit_area
            )

            chamber_pressure_list.append(current_chamber_pressure)

        return chamber_pressure_list

    def evaluate_empirical_burn_rate(self):
        """Computation of the total grain burn rate from user supplied thrust
        data.

        Source: https://pt.scribd.com/document/549653934/RocketElementsHandout.
        Fundamentos de propulsao solida de foguetes. Equation (5.26) on page 42,
        derived from Sutton, Rocket Propulsion Elements, 8ed, page 445 - 453.

        Returns:
            list: burn rate values from thrust data
        """
        self.empirical_chamber_pressure = self.evaluate_empirical_chamber_pressure()

        burn_rate_list = []

        # Initial conditions
        free_volume = self.motor.free_volume
        regressed_length = 0.0
        burn_area = self.motor.total_burn_area

        # Loop iteration - calculate burn_rate at each step
        for pressure_list_index, chamber_pressure in enumerate(
            self.empirical_chamber_pressure[:-1]
        ):
            delta_time = (
                self.empirical_time_steps[pressure_list_index + 1]
                - self.empirical_time_steps[pressure_list_index]
            )
            delta_pressure = (
                self.empirical_chamber_pressure[pressure_list_index + 1]
                - self.empirical_chamber_pressure[pressure_list_index]
            )
            chamber_pressure_derivative = delta_pressure / delta_time

            burn_rate = self.evaluate_burn_rate(
                chamber_pressure,
                chamber_pressure_derivative,
                free_volume,
                burn_area,
            )

            burn_rate_list.append(burn_rate)
            self.grain.regressed_length += burn_rate * delta_time

            if (
                self.grain.regressed_length
                >= self.grain.outer_radius - self.grain.initial_inner_radius
            ):
                break
            burn_area = self.motor.total_burn_area
            free_volume += burn_area * burn_rate * delta_time

        return burn_rate_list


class EmpiricalExport(Export):
    def __init__(self, BurnEmpirical):
        self.BurnEmpirical = BurnEmpirical
        self.post_processing()
        self.data_regression()

    def post_processing(self):
        (
            self.max_empirical_chamber_pressure,
            self.max_empirical_thrust,
        ) = Export.evaluate_max_variables_list(
            self.BurnEmpirical.empirical_time_steps,
            [
                self.BurnEmpirical.empirical_chamber_pressure,
                self.BurnEmpirical.empirical_thrust,
            ],
        )

        self.empirical_total_impulse = self.BurnEmpirical.evaluate_total_impulse(
            self.BurnEmpirical.empirical_thrust,
            self.BurnEmpirical.empirical_time_steps,
        )
        self.empirical_specific_impulse = self.BurnEmpirical.evaluate_specific_impulse(
            self.BurnEmpirical.empirical_thrust,
            self.BurnEmpirical.empirical_time_steps,
        )

        return None

    def data_regression(self):
        """Evaluation of polynomial regression from empirical burnrate.

        Returns:
            None
        """
        self.max_pressure_index = np.argmax(
            self.BurnEmpirical.empirical_chamber_pressure
        )

        self.pressurization_regression = np.polynomial.polynomial.Polynomial.fit(
            self.BurnEmpirical.empirical_chamber_pressure[: self.max_pressure_index],
            self.BurnEmpirical.empirical_burn_rate[: self.max_pressure_index],
            1,
        )

        self.depressurization_regression = np.polynomial.polynomial.Polynomial.fit(
            self.BurnEmpirical.empirical_chamber_pressure[self.max_pressure_index :][
                :-1
            ],
            self.BurnEmpirical.empirical_burn_rate[self.max_pressure_index :],
            1,
        )

        return None

    def all_info(self):
        print("Total Impulse: {:.2f} Ns".format(self.empirical_total_impulse))
        print("Max Thrust: {:.2f} N at {:.2f} s".format(*self.max_empirical_thrust))
        print(
            "Mean Thrust: {:.2f} N".format(np.mean(self.BurnEmpirical.empirical_thrust))
        )
        print(
            "Max Chamber Pressure: {:.2f} bar at {:.2f} s".format(
                self.max_empirical_chamber_pressure[0] / 1e5,
                self.max_empirical_chamber_pressure[1],
            )
        )
        print(
            "Mean Chamber Pressure: {:.2f} bar".format(
                np.mean(self.max_empirical_chamber_pressure) / 1e5
            )
        )
        print("Specific Impulse: {:.2f} s".format(self.empirical_specific_impulse))
        print(
            "Burnout Time: {:.2f} s".format(self.BurnEmpirical.empirical_time_steps[-1])
        )
        return None

    def plotting(self):
        try:
            plt.figure(101, figsize=(16, 9))
            plt.plot(
                self.BurnEmpirical.empirical_time_steps,
                self.BurnEmpirical.empirical_thrust,
                color="g",
                linewidth=0.75,
                label=r"$F^{emp}_T$",
            )
            plt.grid(True)
            plt.xlabel("time (s)")
            plt.ylabel("thrust (N)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Thrust as function of time")
            plt.show()
            plt.close()

            plt.figure(102, figsize=(16, 9))
            plt.plot(
                self.BurnEmpirical.empirical_time_steps,
                self.BurnEmpirical.empirical_chamber_pressure,
                color="g",
                linewidth=0.75,
                label=r"$p^{emp}_c$",
            )
            plt.grid(True)
            plt.xlabel("time (s)")
            plt.ylabel("chamber pressure (pa)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Chamber Pressure as function of time")
            plt.show()
            plt.close()

            plt.figure(103, figsize=(16, 9))
            plt.plot(
                self.BurnEmpirical.empirical_time_steps[:-1],
                self.BurnEmpirical.empirical_burn_rate,
                color="g",
                linewidth=0.75,
                label=r"$r^{emp}$",
            )
            plt.grid(True)
            plt.xlabel("time (s)")
            plt.ylabel("burn rate (m/s)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Burn Rate as function of time")
            plt.show()
            plt.close()

            plt.figure(104, figsize=(16, 9))
            plt.plot(
                *self.pressurization_regression.linspace(100),
                color="g",
                linewidth=0.75,
                label=r"$r^{press}$",
            )
            plt.plot(
                *self.depressurization_regression.linspace(100),
                color="m",
                linewidth=0.75,
                label=r"$r^{depress}$",
            )
            plt.plot(
                self.BurnEmpirical.empirical_chamber_pressure[
                    : self.max_pressure_index
                ],
                self.BurnEmpirical.empirical_burn_rate[: self.max_pressure_index],
                color="black",
                linewidth=0.75,
                linestyle="-.",
            )
            plt.plot(
                self.BurnEmpirical.empirical_chamber_pressure[
                    self.max_pressure_index :
                ][:-1],
                self.BurnEmpirical.empirical_burn_rate[self.max_pressure_index :],
                color="black",
                linewidth=0.75,
                linestyle="-.",
            )
            plt.grid(True)
            plt.xlabel("chamber pressure (pa)")
            plt.ylabel("burn rate (m/s)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Burn Rate as function of chamber pressure")
            plt.show()
            plt.close()

        except AttributeError:
            print(">>> Empirical data not found, empirical plots were not updated.\n")

        return None


if __name__ == "__main__":
    """Burn definitions"""
    Grao_Leviata = Bates(
        outer_radius=71.92 / 2000,
        inner_radius=31.92 / 2000,
    )
    Star_Test = Star(
        outer_radius=71.92 / 2000,
        star_maximum=(71.92 / 2000) / 3 * (5 / 3),
        star_minimum=(71.92 / 2000) / 9 * 2,
        star_points=4,
    )
    Leviata = Motor(
        Star_Test,
        grain_number=4,
        chamber_inner_radius=77.92 / 2000,
        nozzle_throat_radius=17.5 / 2000,
        nozzle_exit_radius=44.44 / 2000,
        nozzle_angle=15 * np.pi / 180,
        chamber_length=600 / 1000,
    )
    KNSB = Propellant(
        specific_heat_ratio=1.1361,
        density=1700,
        products_molecular_mass=39.9e-3,
        combustion_temperature=1600,
        # burn_rate_a=5.13,
        # burn_rate_n=0.22,
        interpolation_list="data/burnrate/KNSB3.csv",
    )

    Ambient = Environment(latitude=-0.38390456, altitude=627, ellipsoidal_model=True)

    """Static-fire data"""

    data_path = "data/static_fires/leviata_final_curve.csv"
    ext_data = np.loadtxt(
        data_path,
        delimiter=",",
        unpack=True,
        skiprows=1,
    )

    Empirical_Simulation = BurnEmpirical(
        Leviata, KNSB, environment=Ambient, empirical_data=ext_data
    )
    ExportPlot = EmpiricalExport(Empirical_Simulation)

    ExportPlot.all_info()
    ExportPlot.plotting()
