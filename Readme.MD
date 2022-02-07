# SolidPy

SolidPy is _Projeto Jupiter_ solid motor simulation code. This repository is under development and is a work in progress.

## Goals

This project aims to build an easy-to-use and versatile simulation tool for solid motor design and validation inside the project. It will follow an OOP design which allows for straightforward simulations of different combinations of propellants, motors and external conditions.

### Expected Outputs

- Total Impulse
- Specific Impulse
- Thrust (time)
- Pressure Chamber (time)
- Kn (time)

## Prerequisites

- Python 3
- NumPy >= 1.0
- SciPy >= 1.0
- Matplotlib >= 3.0

## Authors

-
-

## Assumptions

- Unidirectional and isentropic flow
- Homogeneous combustion products
- No heat exchange with chamber walls (adiabatic)
- No shockwaves or discontinuities in nozzle
- Erosive burning is neglected
- BATES grain
- some others
