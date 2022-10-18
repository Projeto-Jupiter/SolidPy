from setuptools import setup

setup(
    name="solidpy",
    packages=["solidpy"],
    install_requires=[
        "ambiance>=1.3.0",
        "matplotlib>=3.5.1",
        "numpy>=1.19.5",
        "scikit_fmm>=2022.8.15",
        "scipy>=1.7.0"
    ],
)