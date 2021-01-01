"""Standard library imports"""
import os
import numpy as np  # 1.19.4
import pandas as pd  # 1.2.0
import scipy as sp  # 1.6.0

print(f'numpy vers {np.__version__} , \t pandas vers {pd.__version__} , \t scipy vers {sp.__version__}')


# print(os.listdir(os.getcwd()))


class Rotor():
    """Rotor object can be used to calculate the normal and tangencial loads for
    the DTU 10 MW Wind Turbine."""
    # class attributes

    number_of_blades = 3
    radio = 89.17

    """Load the aerodynamics profiles."""

    ffa_241 = pd.read_csv('FFA-W3-241.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_301 = pd.read_csv('FFA-W3-241.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_360 = pd.read_csv('FFA-W3-241.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_480 = pd.read_csv('FFA-W3-241.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_600 = pd.read_csv('FFA-W3-241.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    cylinder = pd.read_csv('FFA-W3-241.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])

    ffa_dict = dict(zip(np.arange(6), [ffa_241, ffa_301, ffa_360, ffa_480, ffa_600, cylinder]))

    """Load the balde data: twist, chord and thick/chord relation for corresponding
    radial position."""

    blade_data = pd.read_csv('bladedat.txt', sep='\t', names=['r', 'twist', 'c', 't/c'])

    # class constructor

    def __init__(self, rho=1.225):
        self.rho = rho # [kg/m3]

    @classmethod
    def get_blade_data(cls):
        return cls.blade_data

    def load_calc(self, wind_speed):
        print(len(Rotor.blade_data))
        # assert isinstance(wind_speed, object)
        v0 = wind_speed
        print(self.rho)
        for row in range(len(Rotor.blade_data)):
            #a, a_prime, k, i, conv, ac, diff = 0, 0, 0, 0, 1, 1 / 3, 0.0001
            r, c, t, twist = Rotor.blade_data['r'][row], Rotor.blade_data['c'][row], Rotor.blade_data['t/c'][row], Rotor.blade_data['twist'][row]
            print(r,c,t,twist)

        return


rotor = Rotor()
rotor.load_calc(wind_speed=9)
