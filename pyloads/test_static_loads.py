import pandas as pd
import numpy as np

from pyloads.static_loads import Rotor
import unittest


class TestStaticLoads(unittest.TestCase):
    def test_normal_tangential_loads(self):
        """Test for normal_tangential_loads results for 6 m/s wind speed."""

        dtu_10mw = Rotor()
        u, pitch, rpm = dtu_10mw.show_operation(u=6)
        tsr = (rpm * np.pi / 30) * Rotor.radio / u

        tan_i, norm_i = dtu_10mw.normal_tangential_loads(tsr, u, dtu_10mw.radio[0], dtu_10mw.twist[0] + pitch,
                                                         dtu_10mw.cord[0], dtu_10mw.t_c[0], verbose=False)

        self.assertEqual(np.round(tan_i, 3), -16.636, 'should be -16.636 ')  # -16.635578670240793
        self.assertEqual(np.round(norm_i, 3), 56.735, 'should be 56.735 ')  # 56.735026640634054

    def test_power_thrust_coefficient(self):
        dtu_10mw = Rotor()
        u, pitch, rpm = dtu_10mw.show_operation(u=6)
        tsr = (rpm * np.pi / 30) * Rotor.radio / u

        power, thrust, pT, pN = dtu_10mw.power_thrust_coefficient(tsr, u, dtu_10mw.radio, dtu_10mw.twist + pitch,
                                                                  dtu_10mw.cord,
                                                                  dtu_10mw.t_c, plot_Loads=False)

        self.assertEqual(np.round(power), 1533727.0, 'should be 1533727.0')
        self.assertEqual(np.round(thrust), 481016.0, 'should be 481016.0')


if __name__ == '__main__':
    unittest.main()

    rot = Rotor()
    Cl, Cd = rot.lift_drag_coeff(alpha=0, t_c=24.1)
    print(f'Cl {Cl}')
    print(f'Cd {Cd}')

    assert rot.lift_drag_coeff(alpha=0, t_c=24.1)[0] == 0.3391, 'review lift_drag_coeff method'
