"""Standard library imports"""
import numpy as np  # 1.19.4
import pandas as pd  # 1.2.0
import scipy as sp  #
from scipy.interpolate import interp1d

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

    blade_data = pd.read_csv('bladedat.txt', sep='\t', names=['r', 'twist', 'c', 't/c'])

    # class constructor

    def __init__(self, rho=1.225):
        self.rho = rho  # [kg/m3]

    @classmethod
    def get_blade_data(cls):
        """Display the blade data: twist, chord and thick/chord relation for corresponding
        radial position."""
        return cls.blade_data

    def lift_drag_coeff(self, t, alpha):
        # TODO:
        #   raise exceptions for thick and alpha
        """Interpolation for drag and lift coefficients.

         Returns: (Cl, Cd)
         :type t: float"""
        cdthick, clthick = np.zeros(6), np.zeros(6)

        for k in range(6):
            f1cl = interp1d(Rotor.ffa_dict[k].iloc[:, 0], Rotor.ffa_dict[k].iloc[:, 1])
            f1cd = interp1d(Rotor.ffa_dict[k].iloc[:, 0], Rotor.ffa_dict[k].iloc[:, 2])
            clthick[k] = f1cl(alpha * 180 / np.pi)  # must convert into degrees
            cdthick[k] = f1cd(alpha * 180 / np.pi)

        thick_prof = np.array([24.1, 30.1, 36., 48., 60., 100.])
        f2cl = interp1d(thick_prof, clthick)
        f2cd = interp1d(thick_prof, cdthick)
        Cl = f2cl(t)
        Cd = f2cd(t)
        print(f'Cd= {Cd} \t Cl= {Cl}')
        return Cl, Cd

    def load_calc(self, wind_speed, pitch, rpm):
        v0 = wind_speed

        try:
            v0 + 2
        except TypeError:
            print('The wind speed should be a number (int or float)')

        assert (v0 > 0), 'The wind speed should be positive'

        print(len(Rotor.blade_data))
        # assert isinstance(wind_speed, object)

        print(self.rho)
        for row in range(len(Rotor.blade_data)):
            # FIXME: complete loads calculation code.
            a, a_prime, k, i, conv, ac, diff = 0, 0, 0, 0, 1, 1 / 3, 0.0001
            r, c, t, twist = Rotor.blade_data['r'][row], Rotor.blade_data['c'][row], Rotor.blade_data['t/c'][row], \
                             Rotor.blade_data['twist'][row]
            while k == 0:
                i += 1
                flow_angle = np.arctan(((1 - a) * v0) / ((1 + a_prime) * w * r))
                alpha = flow_angle - (twist * np.pi / 180 + pitch * np.pi / 180)  # alpha in radians

                Cn = Cl * np.cos(flow_angle) + Cd * np.sin(flow_angle)
                Ct = Cl * np.sin(flow_angle) - Cd * np.cos(flow_angle)

                v_rel = ((1 - a) * v0) / np.sin(flow_angle)

                pN_load = 0.5 * rho * v_rel ** 2 * c * Cn
                pT_load = 0.5 * rho * v_rel ** 2 * c * Ct

                F = (2 / np.pi) * np.arccos(np.exp(-(B / 2) * ((R - r) / (r * np.sin(abs(flow_angle))))))

                sigma = c * B / (2 * np.pi * r)

        return a

    def flow_angle(self):
        # TODO: define flow angle method.
        pass


if __name__ == "__main__":
    rotor = Rotor()
    rotor.load_calc(wind_speed=9)
