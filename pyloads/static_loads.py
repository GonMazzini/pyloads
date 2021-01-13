"""Standard library imports"""
import numpy as np  # 1.19.4
import pandas as pd  # 1.2.0
import scipy as sp  #
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

print(f'numpy version {np.__version__} , \t pandas vers {pd.__version__} , \t scipy vers {sp.__version__}')


# TODO
#   LAST ERROR: ValueError: A value in x_new is above the interpolation range.
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
    ffa_301 = pd.read_csv('FFA-W3-301.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_360 = pd.read_csv('FFA-W3-360.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_480 = pd.read_csv('FFA-W3-480.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    ffa_600 = pd.read_csv('FFA-W3-600.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])
    cylinder = pd.read_csv('cylinder.txt', sep='\t', names=['alpha', 'Cl', 'Cd', 'Cm'])

    ffa_dict = dict(zip(np.arange(6), [ffa_241, ffa_301, ffa_360, ffa_480, ffa_600, cylinder]))

    blade_data = pd.read_csv('bladedat.txt', sep='\t', names=['r', 'twist', 'c', 't/c'])

    # class constructor

    def __init__(self, rho=1.225):
        self.rho = rho  # [kg/m3]

    @classmethod
    def get_blade_data(cls, property='all'):
        """Display the blade data: twist, chord and thick/chord relation for corresponding
        radial position."""

        if property == 'all':
            output = cls.blade_data
        elif property == 'r':
            output = cls.blade_data['r']
        elif property == 'twist':
            output = cls.blade_data['twist']
        elif property == 'c':
            output = cls.blade_data['c']
        elif property == 't/c':
            output = cls.blade_data['t/c']

        return output

    def lift_drag_coeff(self, alpha, t_c):
        # TODO:
        #   raise exceptions for thick and alpha
        #   read how to write docstring in Python PEP-8
        """Interpolation for drag and lift coefficients.

         Returns: (Cl, Cd)

         --------------parameters:

         t_c: float (ie. 24.1)
         alpha: int or float (rad)
         """

        t = t_c
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

        #print(f'alpha:{alpha}, t_c={t_c}')
        #print(f'Cd= {Cd} \t Cl= {Cl}')

        return Cl, Cd

    def normal_tangential_loads(self, tsr, v_0, theta, r, c, t_c, a=0.2, aa=0.2, i=0, imax=100):

        def glauert_equation(x, sigma, F, phi, Cn):
            return [x[0] - ((1 - x[1]) ** 2 * sigma * Cn) / (np.sin(phi) ** 2),
                    x[0] - 4 * x[1] * (1 - 0.25 * (5 - 3 * x[1]) * x[1]) * F]

        tol_a, tol_aa = 10, 10
        B = Rotor.number_of_blades
        sigma = (c * B) / (2 * np.pi * r)

        while tol_a > 10 ** (-3) and tol_aa > 10 ** (-3) and i < imax:
            a0, aa0 = a, aa
            phi = np.arctan(((1 - a) * Rotor.radio) / ((1 + aa) * tsr * r))
            alpha = np.rad2deg(phi) - theta
            alpha = np.deg2rad(alpha)
            Cl, Cd = self.lift_drag_coeff(alpha, t_c)
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
            F = (2 / np.pi) * np.arccos(np.exp(-(B / 2) * (Rotor.radio - r) / (r * np.sin(abs(phi)))))

            if a <= 1 / 3:
                a = 1 / (((4 * F * np.sin(phi) ** 2) / (sigma * Cn)) + 1)

            else:
                CT = fsolve(glauert_equation, [1, a], args=(sigma, F, phi, Cn))  # [1, a] is necessary. why?

            aa = 1 / (((4 * F * np.sin(phi) * np.cos(phi)) / (sigma * Ct)) - 1)
            tol_a, tol_aa = abs(a - a0), abs(aa - aa0)
            i += 1

        v_rel = (v_0 / np.sin(phi)) * (1 - a)
        pT = 0.5 * Ct * self.rho * (v_rel ** 2) * c
        pN = 0.5 * Cn * self.rho * (v_rel ** 2) * c

        if i == imax:
            print('warning: Not converged')

        return pT, pN

    def power(self, tsr, u, theta, r, c, t_c, plot_Loads = False):
        pT = np.zeros(len(r))
        pN = np.zeros(len(r))
        for i in range(len(r)):
            try:
                if i != 0 and i != 1:
                    pass
                else:
                    print(f'tsr={tsr}. u={u}, theta={theta[i]}, r={r[i]}, c={c[i]}, t_c={t_c[i]}')
                pT[i], pN[i] = self.normal_tangential_loads(tsr, u, theta[i], r[i], c[i], t_c[i])
            except TypeError:
                pT[i], pN[i] = np.nan, np.nan
        # append and assign values at r=R
        r = np.append(r, Rotor.radio)
        pT = np.append(pT[:-1], 0)  # The -1 is a rusty way to solve the problem
        pN = np.append(pN[:-1], 0)
        w = tsr * u / Rotor.radio
        # P = integrate(pT, r) * B * w
        # T = thruster(pN, r)


        if plot_Loads:  # ( == True)
            # TODO
            #   Add labels to grahps.
            plt.figure()
            plt.plot(Rotor.get_blade_data(property='r'), pN)
            plt.plot(Rotor.get_blade_data(property='r'), pT)
            plt.grid()
            plt.ylabel('normal loads [N/m]', fontsize=14)
            plt.xlabel('rotor radius [m]', fontsize=14)
            plt.show()

        return pT, pN


if __name__ == "__main__":
    # instance a rotor object.
    WT_data = pd.read_csv('operation.txt', sep='\s+')
    WT_data.index = WT_data.u
    print(WT_data.loc[6])
    u, pitch, rpm = WT_data.loc[6]
    rotor = Rotor()
    print(type(rotor))  # <class '__main__.Rotor'>
    # test power method
    tsr = (rpm * np.pi / 30) * Rotor.radio / u

    A, B = rotor.power(tsr, u, Rotor.blade_data['twist'] + pitch, Rotor.blade_data['r'], Rotor.blade_data['c'],
                Rotor.blade_data['t/c'], plot_Loads=True)

    a, b = rotor.normal_tangential_loads(tsr, u, Rotor.blade_data['twist'][0] + pitch, Rotor.blade_data['r'][0], Rotor.blade_data['c'][0],
                Rotor.blade_data['t/c'][0])

# TODO
#   * improve plots
#   * add power and thrust calculation


print('Finish ran static loads.')
