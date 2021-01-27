"""Standard library imports"""
import numpy as np  # 1.19.4
import pandas as pd  # 1.2.0
import scipy as sp  #
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

"""Local modules"""
from pyloads.blade_data import BladeFeatures
from pyloads.aerodynamic_profiles import AeroProfiles
from pyloads.operation_dtu10mw import Operation

print(f'numpy version {np.__version__} , \t pandas vers {pd.__version__} , \t scipy vers {sp.__version__}')


# TODO
#   LAST ERROR: ValueError: A value in x_new is above the interpolation range.
# print(os.listdir(os.getcwd()))


class Rotor(Operation):
    """Rotor object can be used to calculate the normal and tangential loads for
    the DTU 10 MW Wind Turbine.

    Parameters
    ----------
    radio : array , default=None
    twist : array , default=None
    cord : array , default=None
    t_c : array , default=None
    profiles : array , default='Default'"""
    # class attributes
    number_of_blades = 3
    radio = 89.17  # [m]
    rho = 1.225  # [km/m3]

    # blade_data = pd.read_csv('bladedat.txt', sep='\t', names=['r', 'twist', 'c', 't/c'])

    # class constructor

    def __init__(self, radio=None, twist=None, cord=None, t_c=None, profiles='Default'):
        super().__init__()
        bld = BladeFeatures()
        # Take the parameters from BladeFeatures corresponding to DTU 10MW
        if radio is None:
            self.radio = bld.radio

        if twist is None:
            self.twist = bld.twist

        if cord is None:
            self.cord = bld.cord

        if t_c is None:
            self.t_c = bld.t_c

        """Load the aerodynamics profiles."""
        # TODO:
        #   This should be in a child class from Rotor called DTU-10MW
        if profiles == 'Default':
            aero_prof = AeroProfiles()

            ffa_241 = aero_prof.ffa_241
            ffa_301 = aero_prof.ffa_301
            ffa_360 = aero_prof.ffa_360
            ffa_480 = aero_prof.ffa_480
            ffa_600 = aero_prof.ffa_600
            cylinder = aero_prof.cylinder

            self.ffa_dict = dict(zip(np.arange(6), [ffa_241, ffa_301, ffa_360, ffa_480, ffa_600, cylinder]))

        # TODO:
        #   check that the len of the features is the same
        #   check that the data has valid ranges..

    @staticmethod
    def integrate(y, r):
        """Useful function for numerical integration, used for power"""
        M = 0  # dummy assignment before loop
        for k in range(len(y) - 1):
            A_k = (y[k + 1] - y[k]) / (r[k + 1] - r[k])
            B_k = (y[k] * r[k + 1] - y[k + 1] * r[k]) / (r[k + 1] - r[k])
            M += 1 / 3 * A_k * ((r[k + 1]) ** 3 - (r[k]) ** 3) + 0.5 * B_k * ((r[k + 1]) ** 2 - (r[k]) ** 2)
        return M

    @staticmethod
    def thruster(pN, r):
        # [r] m
        T = 0
        B = Rotor.number_of_blades
        for i in range(len(pN) - 1):
            T += (pN[i + 1] + pN[i]) * 0.5 * (r[i + 1] - r[i])
        return T * B

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
            f1cl = interp1d(self.ffa_dict[k].iloc[:, 0], self.ffa_dict[k].iloc[:, 1])
            f1cd = interp1d(self.ffa_dict[k].iloc[:, 0], self.ffa_dict[k].iloc[:, 2])
            clthick[k] = f1cl(alpha * 180 / np.pi)  # must convert into degrees
            cdthick[k] = f1cd(alpha * 180 / np.pi)

        thick_prof = np.array([24.1, 30.1, 36., 48., 60., 100.])
        f2cl = interp1d(thick_prof, clthick)
        f2cd = interp1d(thick_prof, cdthick)
        Cl = f2cl(t)
        Cd = f2cd(t)

        # print(f'alpha:{alpha}, t_c={t_c}')
        # print(f'Cd= {Cd} \t Cl= {Cl}')

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
        pT = 0.5 * Ct * Rotor.rho * (v_rel ** 2) * c
        pN = 0.5 * Cn * Rotor.rho * (v_rel ** 2) * c

        if i == imax:
            print('warning: Not converged')

        return pT, pN

    def power_thrust_coefficient(self, tsr, u, theta, r, c, t_c, plot_Loads=False):
        """Calculate the power and trhust for given operational parameters.
        This method uses the norma_tangencial_loads in order to calculate the loads."""

        pT = np.zeros(len(r))
        pN = np.zeros(len(r))
        for i in range(len(r)):
            try:
                pT[i], pN[i] = self.normal_tangential_loads(tsr, u, theta[i], r[i], c[i], t_c[i])
            except TypeError:
                pT[i], pN[i] = np.nan, np.nan
        # append and assign values at r=R
        r = np.append(r, Rotor.radio)
        pT = np.append(pT[:-1], 0)  # The -1 is a rusty way to solve the problem
        pN = np.append(pN[:-1], 0)
        w = tsr * u / Rotor.radio
        power = Rotor.integrate(pT, r) * Rotor.number_of_blades * w
        thrust = Rotor.thruster(pN, r)

        if plot_Loads:  # ( == True)
            plt.figure()
            plt.plot(self.radio, pN)
            plt.plot(self.radio, pT)
            plt.grid()
            plt.ylabel('normal loads [N/m]', fontsize=14)
            plt.xlabel('rotor radius [m]', fontsize=14)
            plt.show()

        return power, thrust, pT, pN

    def power_curve(self, u_vector, w_vector, pitch_vector, plot_curve=True):
        """Calculate and plot the Power Curve given a vector of wind speed,
        rotational speed (in RPM) and corresponding pitch angle."""

        P, T = np.zeros(len(u_vector)), np.zeros(len(u_vector))
        #  df = Rotor.blade_data.iloc[0:]
        pN = np.zeros([len(u_vector), len(self.radio)])
        pT = np.zeros([len(u_vector), len(self.radio)])
        print(u_vector)
        for j in range(len(u_vector)):
            u = u_vector.values[j]
            w = w_vector.values[j] * np.pi / 30  # convert from RPM to rad/s
            pitch = pitch_vector.values[j]
            TSR = w * Rotor.radio / u
            print(TSR, w, pitch)
            P[j], T[j], pT[j,], pN[j,] = self.power_thrust_coefficient(TSR, u, self.twist + pitch, self.radio,
                                                                       self.cord, self.t_c)
        if plot_curve:
            plt.plot(u_vector, P / 1e6, linestyle='--', marker='o')
            plt.xlabel('Wind speed')
            plt.ylabel('power [MW]')
            plt.grid()
        return P, T


if __name__ == "__main__":
    # instance a rotor object.
    # WT_data = pd.read_csv('operation.txt', sep='\s+')
    # WT_data.index = WT_data.u
    # print(WT_data.loc[6])
    # u, pitch, rpm = WT_data.loc[6]
    dtu_10mw = Rotor()
    print(type(dtu_10mw))  # <class '__main__.Rotor'>
    # test power method
    u, pitch, rpm = dtu_10mw.show_operation(u=4)
    tsr = (rpm * np.pi / 30) * Rotor.radio / u

    # P, T = rotor.power_curve(WT_data.u, WT_data.RPM, WT_data.pitch)

    power, thrust, pT, pN = dtu_10mw.power_thrust_coefficient(tsr, u, dtu_10mw.twist + pitch, dtu_10mw.radio,
                                                              dtu_10mw.cord,
                                                              dtu_10mw.t_c, plot_Loads=True)

    # tan_i, norm_i = dtu_10mw.normal_tangential_loads(tsr, u, dtu_10mw.twist[0] + pitch,
    #                                              dtu_10mw.radio[0],
    #                                              dtu_10mw.cord[0], dtu_10mw.t_c[0])

# TODO
#   * Review... different results as IPYNB :
#     [56.73502664, 128.66511904, 196.82870201, 955.83255009,
#      1372.0058535, 1685.01982088, 2137.75442846, 2601.22578824,  <<<<< 2601 is the first different value.
#      2986.14621076, 3263.07927127, 3431.68926186, 3461.78866371,
#      3474.49291476, 3421.50774457, 3224.94102283, 2815.05298268,
#      2063.40641495, 0.])
#   * Improve plots
#   * let DTU_10MW be a subclass and define Rotor with user params.
#   * add power and thrust calculation DONE !
#   * add DEFLECTION !


print('Finish ran static loads.')
