# TODO
#   CONTINUE WITH TRANSMAT CODE.
#   CHECK CONSISTENCY WITH MATLAB CODE (DIFF RESULTS)




class TransMat:

    def __init__(self):
        import numpy as np
        from pyloads.dtu10mw.blade_data import BladeFeatures
        bld = BladeFeatures()
        self.r = bld.radio
        self.cord = bld.cord
        self.twist = bld.twist
        self.t_c = bld.t_c


        self.vo_hub = 10  # wind speed at hub height

        self.omega = 0.62  # rotational velocity [rad/s]

        self.theta_cone = 0.  # angle between rotor plane and blade [rad]
        self.theta_yaw = 0. * np.pi / 180  # angle between incoming wind speed and shaft
        self.theta_tilt = 0.  # elevation angle of the nacelle
        self.theta_b1 = np.zeros(360)

        self.diamater = 178.28

        self.H = 119
        self.L = 7.1

        self.bladepos = 70

        self.time_rot = 2 * np.pi / self.omega  # seconds to realize one rotation
        self.dtheta = np.pi / 180
        self.dt = self.dtheta / self.omega

        self.nu = 0  # shear exp
        self.a = 2.5  # tower radius

        # class atributes

    def time_independent_mat(self):
        """ Return the time independent matrix a21 (from sys2 to sys1) and matrix a34 (from sys3 to sys 4."""

        a1_yaw = np.array([[1, 0, 0],
                           [0, np.cos(self.theta_yaw), np.sin(self.theta_yaw)],
                           [0, -np.sin(self.theta_yaw), np.cos(self.theta_yaw)]])

        a2_tilt = np.array([[np.cos(self.theta_tilt), 0, -np.sin(self.theta_tilt)],
                            [0, 1, 0],
                            [np.sin(self.theta_tilt), 0, np.cos(self.theta_tilt)]])

        a3_cone = np.array([[1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]])

        a12 = np.dot(np.dot(a1_yaw, a2_tilt), a3_cone)
        a21 = np.transpose(a12)

        # Compute matrix from sys3 to sys4.

        a34 = np.array([[np.cos(self.theta_cone), 0, -np.sin(self.theta_cone)],
                        [0, 1, 0],
                        [np.sin(self.theta_cone), 0, np.cos(self.theta_cone)]])

        return a12, a34

    def time_loop(self):
        time = np.zeros(361)  # will use 360 time stamps
        # pre allocate
        theta = np.zeros([360,3])
        rpos_b1_x, rpos_b1_y, rpos_b1_z = np.zeros(360), np.zeros(360), np.zeros(360)


        # define coordinates in system 1 and system 2.
        r_tow_sys1 = np.array([H, 0, 0])
        r_shaft_sys2 = np.array([0, 0, -L])

        for nt in range(1, 360):
            time[nt] = (nt - 1) * self.dt
            theta[nt, 0] = theta[nt - 1,0] + self.dtheta
            theta[nt, 1] = theta[nt, 0] + 2 * np.pi / 3
            theta[nt, 2] = theta[nt, 0] + 4 * np.pi / 3

            # calculate matrix a23 for each blade. Note these are TIME DEPENDENT.
            for b in range(3):
                a23 = np.array([[np.cos(theta[nt, b]), np.sin(theta[nt, b]), 0],
                                            [-np.sin(theta[nt, b]), np.cos(theta[nt, b]), 0],
                                            [0, 0, 1]])

                # transformation matrix from system 1 to 4
                # call the time_independent method to get a12 and a34.
                a12, a34 = self.time_independent_mat()

                # Blade 1
                a14 = np.dot(a34, np.dot(a23, a12))
                ab41 = np.transpose(a14)

                # loop for each blade position
                for i in range(len(self.r)-1):



                



            # Position in original system
            r_tow_sys1 = np.array([self.H, 0, 0])
            r_shaft_sys2 = np.array([0, 0, -1])
            r_blade_sys4 = np.array([self.bladepos, 0, 0])

            # Positions in blade system for blade 1
            #rpos_b1_sys1 = r_tow_sys1 + a21 * r_shaft_sys2 + a41_B1 * r_blade_sys4

        print("stop")
        return True


if __name__ == "__main__":
    import numpy as np

    wt = TransMat()
    # a21, a34 = wt.time_independent_mat()
    # a21, a34 = wt.time_independent_mat()
    tl = wt.time_loop()
