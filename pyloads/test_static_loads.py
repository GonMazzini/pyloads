import pandas as pd
import numpy as np

from pyloads.static_loads import Rotor
import unittest


class TestStaticLoads(unittest.TestCase):
    def test_normal_tangential_loads(self):
        dtu_10mw = Rotor()

        WT_data = pd.read_csv('operation.txt', sep='\s+')
        WT_data.index = WT_data.u
        u, pitch, rpm = WT_data.loc[6]
        tsr = (rpm * np.pi / 30) * Rotor.radio / u

        tan_i, norm_i = dtu_10mw.normal_tangential_loads(tsr, u, dtu_10mw.twist[0] + pitch,
                                                         dtu_10mw.radio[0],
                                                         dtu_10mw.cord[0], dtu_10mw.t_c[0])
        self.assertEqual(tan_i, -16.635578670240793, 'should be -16.635578670240793 ')  # -16.635578670240793
        self.assertEqual(norm_i, 56.735026640634054, 'should be 56.735026640634054 ')  # 56.735026640634054

    def test_power_thrust_coefficient(self):
        dtu_10mw = Rotor()
        print(dtu_10mw.radio)

        WT_data = pd.read_csv('operation.txt', sep='\s+')
        WT_data.index = WT_data.u
        u, pitch, rpm = WT_data.loc[6]
        tsr = (rpm * np.pi / 30) * Rotor.radio / u

        power, thrust, pT, pN = dtu_10mw.power_thrust_coefficient(tsr, u, dtu_10mw.radio, dtu_10mw.twist + pitch,
                                                                  dtu_10mw.cord,
                                                                  dtu_10mw.t_c, plot_Loads=True)

        self.assertEqual(power, 1599385.8391734336, 'should be 1599385.8391734336')
        self.assertEqual(thrust, 501337.71234120766, 'should be 501337.71234120766')


if __name__ == '__main__':
    unittest.main()

    rot = Rotor()
    Cl, Cd = rot.lift_drag_coeff(alpha=0, t_c=24.1)
    print(f'Cl {Cl}')
    print(f'Cd {Cd}')

    assert rot.lift_drag_coeff(alpha=0, t_c=24.1)[0] == 0.3391, 'review lift_drag_coeff method'
