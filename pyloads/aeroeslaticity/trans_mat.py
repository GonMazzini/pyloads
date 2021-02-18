# TODO
#   CONTINUE WITH TRANSMAT CODE.

class TransMat:

    def __init__(self):

        import numpy as np

        self.vo_hub = 10  # wind speed at hub height

        self.omega = 0.62  # rotational velocity [rad/s]

        self.theta_cone = 0  # angle between rotor plane and blade [rad]
        self.theta_yaw = 0.1*np.pi/180  # angle between incoming wind speed and shaft
        self.theta_tilt = 0.1  # elevation angle of the nacelle
        self.theta_b1 = np.zeros(360)

        self.diamater = 178.28

        self.H = 119
        self.L = 7.1

        self.bladepos = 70

        self.time_rot = 2* np.pi / self.omega  # seconds to realize one rotation
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

        return a21, a34

    def time_loop(self):
        time = np.zeros(361)
        theta_b1, theta_b2, theta_b3 = np.zeros(360), np.zeros(360),np.zeros(360)
        for nt in range(1, 360):

            time[nt] = (nt-1) * self.dt
            theta_b1[nt] = theta_b1[nt - 1] + self.dtheta
            theta_b2[nt] = theta_b1[nt] + 2*np.pi/3
            theta_b3[nt] = theta_b1[nt] + 4*np.pi/3

            # calculate matrix a23 for each blade. Note these are TIME DEPENDENT.

            a23_b1 = np.array([[np.cos(theta_b1[nt]), np.sin(theta_b1[nt]), 0],
                           [-np.sin(theta_b1[nt]), np.cos(theta_b1[nt]), 0],
                           [0, 0, 1]])
            a23_b2 = np.array([[np.cos(theta_b2[nt]), np.sin(theta_b2[nt]), 0],
                               [-np.sin(theta_b2[nt]), np.cos(theta_b2[nt]), 0],
                               [0, 0, 1]])

            a23_b3 = np.array([[np.cos(theta_b3[nt]), np.sin(theta_b3[nt]), 0],
                               [-np.sin(theta_b3[nt]), np.cos(theta_b3[nt]), 0],
                               [0, 0, 1]])

        print("stop")
        return True

if __name__ == "__main__":
    import numpy as np

    wt = TransMat()
    a21, a34 = wt.time_independent_mat()
    tl= wt.time_loop()