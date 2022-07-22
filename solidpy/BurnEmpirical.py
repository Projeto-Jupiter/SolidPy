# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import math

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties

from Grain import Grain
from Motor import Motor
from Environment import Environment
from Propellant import Propellant
from Burn import Burn
from Export import Export


class BurnEmpirical(Burn):
    def __init__(
        self,
        grain,
        motor,
        propellant,
        environment=Environment(),
        empirical_data=None,
    ):
        Burn.__init__(self, grain, motor, propellant, environment)
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
        if self.grain.density is not None:
            return self.grain.density
        return self.propellant.density

    def evaluate_empirical_chamber_pressure(self):
        """Computation of motors chamber pressure from user supplied
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
        burn_area = (
            self.motor.grain_number
            * self.motor.grain.evaluate_tubular_burn_area(regressed_length)
        )

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
            regressed_length += burn_rate * delta_time
            burn_area = (
                self.motor.grain_number
                * self.motor.grain.evaluate_tubular_burn_area(regressed_length)
            )
            free_volume += burn_area * regressed_length

        return burn_rate_list


class EmpiricalExport(Export):
    def __init__(self, BurnEmpirical):
        self.BurnEmpirical = BurnEmpirical
        self.post_processing()

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
            plt.figure(1, figsize=(16, 9))
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
            plt.savefig("data/burn_simulation/graphs/empirical_thrust.png", dpi=200)

            plt.figure(2, figsize=(16, 9))
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
            plt.savefig(
                "data/burn_simulation/graphs/empirical_chamber_pressure.png", dpi=200
            )

            k = np.polynomial.polynomial.Polynomial.fit(self.BurnEmpirical.empirical_time_steps[:-1], self.BurnEmpirical.empirical_burn_rate, 7)
            k1 = []
            t1 = []
            for l in self.BurnEmpirical.empirical_time_steps[:-1]:
                if k(l) >= 0:
                    k1.append(k(l))
                    t1.append(l)
                else:
                    break

            plt.figure(3, figsize=(16, 9))
            plt.plot(
                t1,
                k1,
                color="g",
                linewidth=0.75,
                label=r"$r^{emp}$",
            )
            plt.grid(True)
            plt.xlabel("time (s)")
            plt.ylabel("burn rate (m/s)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Burn Rate as function of time")
            plt.savefig("data/burn_simulation/graphs/empirical_burn_rate.png", dpi=200)


            max_time_index=0
            for count, value in enumerate(self.BurnEmpirical.empirical_chamber_pressure):
                if value > 3.5e6:
                    max_time_index = count
                    break

            bur = np.polynomial.polynomial.Polynomial.fit(self.BurnEmpirical.empirical_chamber_pressure[:max_time_index], self.BurnEmpirical.empirical_burn_rate[:max_time_index], 3)

            bur_vec = []
            pvec = [i for i in range(int(0.1e5), int(4.2e6), int(1e4))]

            for i in pvec:
                bur_vec.append(bur(i))

            plt.figure(4, figsize=(16, 9))
            plt.plot(
                self.BurnEmpirical.empirical_chamber_pressure[:max_time_index],
                self.BurnEmpirical.empirical_burn_rate[:max_time_index],
                color="g",
                linewidth=0.75,
                label=r"$r^{emp}$",
            )
            plt.grid(True)
            plt.xlabel("chamber pressure (pa)")
            plt.ylabel("burn rate (m/s)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Burn Rate as function of chamber pressure")
            plt.savefig("data/burn_simulation/graphs/empirical_burn_rate_to_pressure.png", dpi=200)

            plt.figure(5, figsize=(16, 9))
            plt.plot(
                pvec,
                bur_vec,
                color="g",
                linewidth=0.75,
                label=r"$r^{emp}$",
            )
            plt.grid(True)
            plt.xlabel("chamber pressure (pa)")
            plt.ylabel("burn rate (m/s)")
            plt.legend(prop=FontProperties(size=16))
            plt.title("Empirical Burn Rate as function of chamber pressure")
            plt.savefig("data/burn_simulation/graphs/empirical_burn_rate_to_pressure2.png", dpi=200)

            p1 = []
            tt = []
            for i in self.BurnEmpirical.empirical_chamber_pressure[:max_time_index]: p1.append(i*1e-6)
            for t in self.BurnEmpirical.empirical_burn_rate[:max_time_index]: tt.append(t*1e3)

            p2 = []
            tt2 = []
            for i in pvec: p2.append(i*1e-6)
            for t in bur_vec: tt2.append(t*1e3)

            Export.raw_simulation_data_export([p1, tt], "data/burnrate/test_emp.csv", ["Chamber Pressure (MPa)", "Burn rate (mm/s)"])
            Export.raw_simulation_data_export([p2, tt2], "data/burnrate/test_emp2.csv", ["Chamber Pressure (MPa)", "Burn rate (mm/s)"])


        except AttributeError:
            print(">>> Empirical data not found, empirical plots were not updated.\n")

        return None


if __name__ == "__main__":
   """Burn definitions"""

   Grao_Mandioca = Grain(
        outer_radius= 94 / 2000, initial_inner_radius= 32 / 2000, 
        #mass=700 / 1000,
        initial_height = 156 / 1000
    )
   Mandioca = Motor(
        Grao_Mandioca,
        grain_number=5,
        chamber_inner_radius=98 / 2000,
        nozzle_throat_radius=11.4 / 1000,
        nozzle_exit_radius=33.5 / 1000,
        nozzle_angle=15 * np.pi / 180,
        chamber_length=840 / 1000,
    )
   KNSB = Propellant(
        specific_heat_ratio=1.1361,
        density=1600,
        products_molecular_mass=39.9e-3,
        combustion_temperature=1600,
        #burn_rate_a=5.13,
        #burn_rate_n=0.22,
        interpolation_list="data/burnrate/KNSB3.csv",
    )
   Ambient = Environment(latitude=-0.38390456, altitude=627, ellipsoidal_model=True)

   """Static-fire data"""

   data_path = "data/static_fires/mandiocaSF.csv"
   ext_data = np.loadtxt(
        data_path,
        delimiter=",",
        unpack=True,
        skiprows=0,
    )

   Empirical_Simulation = BurnEmpirical(
        Grao_Mandioca, Mandioca, KNSB, empirical_data=ext_data
    )
   ExportPlot = EmpiricalExport(Empirical_Simulation)

   ExportPlot.all_info()
   ExportPlot.plotting()
