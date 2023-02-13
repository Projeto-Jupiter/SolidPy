# -*- coding: utf-8 -*-

"""
SolidPy is Projeto Jupiter propulsion's team attempt to create a sophisticated
internal ballistics simulator for solid propellant rocket engines, with versatility
for different configurations and able to generate all the required data for designed motors.
"""

__author__ = (
    "João Lemes Gribel Soares",
    "Pedro Henrique Marinho Bressan",
    "Thais Silva Melo",
)
__copyright__ = "Copyright 20XX, Projeto Jupiter"
__credits__ = (
    "João Lemes Gribel Soares",
    "Pedro Henrique Marinho Bressan",
    "Thais Silva Melo",
)
__license__ = "MIT"
__version__ = ""
__maintainer__ = "João Lemes Gribel Soares"
__email__ = "jgribel@usp.br"
__status__ = "Production"

import csv
import math
import warnings
import time
import skfmm

import numpy as np
from scipy import integrate
from scipy import linalg
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import scipy.constants as const
import matplotlib.pyplot as plt

from .Export import Export
from .Grain import Bates, Star, CustomGeometry
from .Motor import Motor
from .Propellant import Propellant
from .Burn import Burn, BurnSimulation, BurnExport
from .BurnEmpirical import BurnEmpirical, EmpiricalExport
from .Environment import Environment
from .Rail import Rail, RailExport
