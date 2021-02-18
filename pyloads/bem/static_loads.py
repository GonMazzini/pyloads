"""Standard library imports"""
import numpy as np  # 1.19.4
import pandas as pd  # 1.2.0
import scipy as sp  #
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve

"""Local modules"""
from pyloads.dtu10mw.blade_data import BladeFeatures
from pyloads.dtu10mw.aerodynamic_profiles import AeroProfiles
from pyloads.dtu10mw.operation_dtu10mw import Operation



class Rotor(Operation):
    """Create a Rotor object (by Default is the DTU 10 MW) by defining the blade geometry (cord, twist and
    thickness), as well as the aerodynamic profiles (i.e as output from XFOIL interactive program for design).
    - Calculate the Normal and Tangential loads for given operational conditions.
    - Calculate the Power and Thrust coefficient.
    - Calculate and plot the Rotor Power Curve.

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

    operation = pd.DataFrame(
        {'u': [4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
               22., 23., 24., 25., ],
         'pitch': [2.751, 1.966, 0.896, 0., 0., 0., 0., 0., 4.502, 7.266,
                   9.292, 10.958, 12.499, 13.896, 15.2, 16.432, 17.618, 18.758, 19.86,
                   20.927,
                   21.963, 22.975],
         'RPM': [6., 6., 6., 6., 6.426, 7.229, 8.032, 8.836, 9.6, 9.6, 9.6, 9.6,
                 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, ]
         })

    # blade_data = pd.read_csv('bladedat.txt', sep='\t', names=['r', 'twist', 'c', 't/c'])
    # class constructor

    def __init__(self, radio='DTU_10MW', twist='DTU_10MW', cord='DTU_10MW', t_c='DTU_10MW', profiles='DTU_10MW'):
        super().__init__()
        bld = BladeFeatures()
        # Take the parameters from BladeFeatures corresponding to DTU 10MW
        if radio == 'DTU_10MW':
            self.radio = bld.radio

        if twist == 'DTU_10MW':
            self.twist = bld.twist

        if cord == 'DTU_10MW':
            self.cord = bld.cord

        if t_c == 'DTU_10MW':
            self.t_c = bld.t_c

        """Load the aerodynamics profiles."""
        # TODO:
        #   This should be in a child class from Rotor called DTU-10MW
        if profiles == 'DTU_10MW':
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
        """Useful function for numerical integration.

        parameters
        ----------
        y : array
            function to be integrated
        r : array
            integrate over r radial positions.
        """
        M = 0  # dummy assignment before loop
        for k in range(len(y) - 1):
            A_k = (y[k + 1] - y[k]) / (r[k + 1] - r[k])
            B_k = (y[k] * r[k + 1] - y[k + 1] * r[k]) / (r[k + 1] - r[k])
            M += 1 / 3 * A_k * ((r[k + 1]) ** 3 - (r[k]) ** 3) + 0.5 * B_k * ((r[k + 1]) ** 2 - (r[k]) ** 2)
        return M

    @staticmethod
    def thruster(pN, r):
        """Compute the total Thrust (T * num_blades) for the Rotor.

        parameter
        ---------
        pN : array
            normal loads vector (i.e as returned by norma_tangential_loads method.
        r : array
            vector with radial position [m]"""
        # [r] m
        T = 0
        B = Rotor.number_of_blades
        for i in range(len(pN) - 1):
            T += (pN[i + 1] + pN[i]) * 0.5 * (r[i + 1] - r[i])
            # print(f'thrust {T} for item num {i}')
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

    def normal_tangential_loads(self, tsr, v_0, r, theta, c, t_c, a=0.2, aa=0.2, imax=100, verbose=False):
        # TODO:
        #   Add velocity triangle vector plot
        """Calculate the Tangential and Normal loads (in N/m) given the wind speed and tip speed ratio (tsr) for a given
        radial position (twist, cord length and thickness/cord ratio at the radial position must be given).

        Parameters
        ----------
        tsr : int or float
            Tip speed ratio.
        v_0 : int or float
            Wind speed [m/s].
        r : int or float
            Radial position.
        theta : int or float
            Local twist at r.
        c : int or float
            Cord length.
        t_c : int or float
            Thickness/cord ratio.
        a : int or float, default=0.2
            tangential induction factor
        aa : int or float, default=0.2
            normal induction factor
        imax : int, default=100
            max number of iter before stoping loop.
        verbose : bool, default=False
            if True, print and plot the variables for each iteration.

        Returns
        -------
        pT : float
            tangential load [N/m] at radial position.
        pN : float
            normal load [N/m] at radial position.


        """

        def glauert_equation(x, sigma, F, phi, Cn):
            return [x[0] - ((1 - x[1]) ** 2 * sigma * Cn) / (np.sin(phi) ** 2),
                    x[0] - 4 * x[1] * (1 - 0.25 * (5 - 3 * x[1]) * x[1]) * F]
        i = 0
        tol_a, tol_aa = 10, 10
        B = Rotor.number_of_blades
        sigma = (c * B) / (2 * np.pi * r)
        if verbose:
            a_list, aa_list, phi_list, alpha_list, i_list = [], [], [], [], []

        while tol_a > 10 ** (-3) and tol_aa > 10 ** (-3) and i < imax:
            a0, aa0 = a, aa
            phi = np.arctan(((1 - a) * Rotor.radio) / ((1 + aa) * tsr * r))
            alpha = np.rad2deg(phi) - theta
            alpha = np.deg2rad(alpha)
            Cl, Cd = self.lift_drag_coeff(alpha, t_c)
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)

            # if i == 0: #get info of first values to check
            #     print(Rotor.radio)
            #     print(r)
            #     print(tsr)
            #     print('phi:',phi)
            #     print('theta:',theta)
            #     print('alpha:',alpha)
            #     print('Cl:',Cl)
            #     print('Cl:', Cd)

            F = (2 / np.pi) * np.arccos(np.exp(-(B / 2) * (Rotor.radio - r) / (r * np.sin(abs(phi)))))

            if a <= 1 / 3:
                a = 1 / (((4 * F * np.sin(phi) ** 2) / (sigma * Cn)) + 1)

            else:
                CT, a = fsolve(glauert_equation, [1, a], args=(sigma, F, phi, Cn))  # [1, a] is necessary. why?

            aa = 1 / (((4 * F * np.sin(phi) * np.cos(phi)) / (sigma * Ct)) - 1)
            tol_a, tol_aa = abs(a - a0), abs(aa - aa0)

            if verbose:
                print('iter #:',i)
                print('\t a: ',a)
                print('\t a_prime: ',aa)
                print('\t phi: ',phi)
                print('\t alpha: ',alpha)

                a_list.append(a)
                aa_list.append(aa)
                phi_list.append(phi)
                alpha_list.append(alpha)
                i_list.append(i)

            i += 1

        if verbose:
            print('final iteration (i):',i)
            if i>1:
            # TODO
            #   review if figsize is correct.
                fig, axes = plt.subplots(2,2, figsize=(10,4))
                axes[0, 0].plot(i_list, a_list, marker='o')
                axes[0, 0].set_ylabel('a', fontsize=14)
                axes[0, 0].set_xlabel('iteration num')
                axes[0, 1].plot(i_list, aa_list,marker='o')
                axes[0, 1].set_ylabel('a\'', fontsize=14)
                axes[0, 1].set_xlabel('iteration num')
                axes[1, 0].plot(i_list, phi_list,marker='o')
                axes[1, 0].set_ylabel('phi', fontsize=14)
                axes[1, 0].set_xlabel('iteration num')
                axes[1, 1].plot(i_list, alpha_list,marker='o')
                axes[1, 1].set_ylabel('alpha', fontsize=14)
                axes[1, 1].set_xlabel('iteration num')
                fig.tight_layout(pad=3.0)
                plt.show()

        v_rel = (v_0 / np.sin(phi)) * (1 - a)
        pT = 0.5 * Ct * Rotor.rho * (v_rel ** 2) * c
        pN = 0.5 * Cn * Rotor.rho * (v_rel ** 2) * c

        if i == imax:
            print(f'warning: Not converged for {imax} iter at radial position = {r} m')

        return pT, pN

    def power_thrust_coefficient(self, tsr, u, r, theta, c, t_c, plot_Loads=False):
        """Calculate the power and thrust for given operational parameters.
        This method uses the norma_tangential_loads in order to calculate the loads.

        :parameter
        ----------
        tsr : int or float
            Tip speed ratio.
        u : int or float
            Wind speed [m/s].
        r : int or float
            Radial position.
        c : int or float
            Cord length.
        t_c : int or float
            Thickness/cord ratio.
        plot_Loads : bool, default=False
            Plot the normal and tangential loads against radial position.

        :return
        -------
        power : float
            Total power considering Rotor.number_of_blades.
        thrust : float
            Total thrust considering Rotor.number_of_blades.
        pT : float
            tangential load [N/m] at radial position.
        pN : float
            normal load [N/m] at radial position.

        """

        pT = np.zeros(len(r))
        pN = np.zeros(len(r))
        for i in range(len(r)):
            try:
                pT[i], pN[i] = self.normal_tangential_loads(tsr, u, r[i], theta[i], c[i], t_c[i])
            except TypeError:
                pT[i], pN[i] = np.nan, np.nan
        # append and assign values at r=R
        r = np.append(r, Rotor.radio)
        pT = np.append(pT[:-1], 0)  # The -1 is a rusty way to solve the problem
        pN = np.append(pN[:-1], 0)
        w = tsr * u / Rotor.radio
        power = Rotor.integrate(pT, r) * Rotor.number_of_blades * w
        # print(f'power integral{Rotor.integrate(pT,r)}')
        # print(f'thrust integral{Rotor.thruster(pN, r)}')
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
        rotational speed (in RPM) and corresponding pitch angle.

        :parameter
        ----------
        u_vector : array
            range of speed for the power curve
        w_vector : array
            vector with rotational speed (in RPM)
        pitch_vector : array
            pitch angle vector

        :return
        -------
        """

        P, T = np.zeros(len(u_vector)), np.zeros(len(u_vector))
        #  df = Rotor.blade_data.iloc[0:]
        pN = np.zeros([len(u_vector), len(self.radio)])
        pT = np.zeros([len(u_vector), len(self.radio)])
        # print(u_vector)
        for j in range(len(u_vector)):
            u = u_vector.values[j]
            w = w_vector.values[j] * np.pi / 30  # convert from RPM to rad/s
            pitch = pitch_vector.values[j]
            TSR = w * Rotor.radio / u
            # print(TSR, w, pitch)
            P[j], T[j], pT[j,], pN[j,] = self.power_thrust_coefficient(TSR, u, self.radio, self.twist + pitch,
                                                                       self.cord, self.t_c)


        if plot_curve:
            plt.plot(u_vector, P / 1e6, linestyle='--', marker='o')
            plt.xlabel('Wind speed')
            plt.ylabel('power [MW]')
            plt.grid()
        return P, T


if __name__ == "__main__":
    print(f'numpy version {np.__version__} , \t pandas vers {pd.__version__} , \t scipy vers {sp.__version__}')
    print('stop')
    # instance a rotor object.
    # WT_data = pd.read_csv('operation.txt', sep='\s+')
    # WT_data.index = WT_data.u
    # print(WT_data.loc[6])
    # u, pitch, rpm = WT_data.loc[6]
    dtu_10mw = Rotor()
    print(type(dtu_10mw))  # <class '__main__.Rotor'>
    # test power method
    oper_df = dtu_10mw.show_operation()  # returns a DataFrame
    u, pitch, rpm = dtu_10mw.show_operation(u=6)
    tsr = (rpm * np.pi / 30) * Rotor.radio / u

    P, T = dtu_10mw.power_curve(oper_df.u, oper_df.RPM, oper_df.pitch, plot_curve=True)

    # power, thrust, pT, pN = dtu_10mw.power_thrust_coefficient(tsr, u, dtu_10mw.radio, dtu_10mw.twist + pitch,
    #                                                          dtu_10mw.cord,
    #                                                          dtu_10mw.t_c, plot_Loads=True)

    # tan_i, norm_i = dtu_10mw.normal_tangential_loads(tsr, u, dtu_10mw.radio[0], dtu_10mw.twist[0] + pitch,
    #                                                  dtu_10mw.cord[0], dtu_10mw.t_c[0], verbose=True)

# TODO#
#  * Power and Thrust are NEGATIVE...
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
