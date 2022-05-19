# -*- coding: utf-8 -*-

import numpy as np

from Grain import Grain
from Motor import Motor
from Propellant import Propellant
from Environment import Environment
from Burn import Burn, BurnSimulation


class TestBurn:

    simulation_data = np.loadtxt(
        "data/burn_simulation/burn_data.csv", delimiter=",", unpack=True, skiprows=1
    )
    (
        time,
        thrust,
        chamber_pressure,
        exit_pressure,
        exit_velocity,
        free_volume,
        regressed_length,
    ) = simulation_data

    """Test output simulation range and lenght"""

    def time_range(self):
        for time_steps in self.time:
            assert time_steps >= 0

    def test_thrust_range(self):
        for thr in self.thrust:
            assert -9.81 * 12 < thr < 1400

    def test_chamber_pressure_range(self):
        for press in self.chamber_pressure:
            assert 0 < press < 40e5

    def test_exit_pressure_range(self):
        for chamber_press, exit_press in zip(self.chamber_pressure, self.exit_pressure):
            assert 0 < exit_press < chamber_press

    def test_exit_velocity_range(self):
        for exit_vel in self.exit_velocity:
            assert exit_vel > 0

    def test_free_volume(self):
        for free_vol in self.free_volume:
            assert 0 < free_vol < Leviata.chamber_volume

    def test_regressed_length(self):
        for regressed_len in self.regressed_length:
            assert (
                regressed_len
                < Grao_Leviata.outer_radius - Grao_Leviata.initial_inner_radius
            )

    def test_output_list_size(self):
        lenght = len(self.time)
        for out_list in self.simulation_data:
            assert len(out_list) == lenght

    """Test burn methods"""

    def test_nozzle_mass_flow(self):
        assert simulation_burn.evaluate_nozzle_mass_flow(0) == 0.0


"""Rocket definitions"""

Grao_Leviata = Grain(
    outer_radius=71.92 / 2000, initial_inner_radius=31.92 / 2000, mass=700 / 1000
)
Leviata = Motor(
    Grao_Leviata,
    grain_number=4,
    chamber_inner_radius=77.92 / 2000,
    nozzle_throat_radius=17.5 / 2000,
    nozzle_exit_radius=44.44 / 2000,
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

Ambient = Environment(101325, 1.25, -0.38390456)

simulation_burn = BurnSimulation(Grao_Leviata, Leviata, KNSB, Ambient)

"""Static-fire data"""
data_path = "data/static_fires/leviata_raw_data-1.csv"
ext_data = np.loadtxt(
    data_path,
    delimiter=",",
    unpack=True,
    skiprows=1,
)
