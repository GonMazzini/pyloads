"""Standard library imports"""
import numpy as np  # 1.19.4
import pandas as pd  # 1.2.0
import scipy as sp  #
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
print(f'numpy version {np.__version__} , \t pandas vers {pd.__version__} , \t scipy vers {sp.__version__}')


# print(os.listdir(os.getcwd()))


class Rotor():
    """Rotor object can be used to calculate the normal and tangencial loads for
    the DTU 10 MW Wind Turbine."""
    # class attributes
    number_of_blades = 3
    radio = 89.17

    """Load the aerodynamics profiles."""
    # TODO:
    #   This should be in a child class from Rotor called DTU-10MW
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
        #   read how to write docstring in Python PEP-8
        """Interpolation for drag and lift coefficients.

         Returns: (Cl, Cd)
         :type t: float (ie. 24.1)
         """
        # :type alpha: float (radian)

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

    def load_calc(self, TSR,v_0,theta,r,c,t_c,a=0.2,aa=0.2,i=0,imax=100):
        tol_a, tol_aa = 10, 10
        sigma = (c * B) / (2 * m.pi * r)

        while tol_a > 10 ** (-3) and tol_aa > 10 ** (-3) and i < imax:
            a0, aa0 = a, aa
            phi = m.atan(((1 - a) * R) / ((1 + aa) * TSR * r))
            alpha = np.rad2deg(phi) - theta
            Cl, Cd = Double_interpol(alpha, t_c)
            Cn = Cl * m.cos(phi) + Cd * m.sin(phi)
            Ct = Cl * m.sin(phi) - Cd * m.cos(phi)
            F = (2 / m.pi) * m.acos(m.exp(-(B / 2) * (R - r) / (r * m.sin(abs(phi)))))

            if a <= 1 / 3:
                a = 1 / (((4 * F * m.sin(phi) ** 2) / (sigma * Cn)) + 1)

            else:
                CT, a = fsolve(Glauert_eq, [1, a], args=(sigma, F, phi, Cn))

            aa = 1 / (((4 * F * m.sin(phi) * m.cos(phi)) / (sigma * Ct)) - 1)
            tol_a, tol_aa = abs(a - a0), abs(aa - aa0)
            i += 1

        v_rel = (v_0 / m.sin(phi)) * (1 - a)

        return a

    def flow_angle(self):
        # TODO: define flow angle method.
        pass


if __name__ == "__main__":
    rotor = Rotor()
    rotor.lift_drag_coeff(t=24.1,alpha=2*np.pi/180)

    print('Finish ran static loads.')