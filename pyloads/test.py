from pyloads.static_loads import Rotor
import unittest


class TestStaticLoads(unittest.TestCase):
    def test_normal_tangential_loads(self):
        dtu_10mw = Rotor()
        import pandas as pd
        import numpy as np
        WT_data = pd.read_csv('operation.txt', sep='\s+')
        WT_data.index = WT_data.u
        u, pitch, rpm = WT_data.loc[6]
        tsr = (rpm * np.pi / 30) * Rotor.radio / u
        tan_i, norm_i = dtu_10mw.normal_tangential_loads(tsr, u, Rotor.blade_data['twist'][0] + pitch,
                                                         Rotor.blade_data['r'][0],
                                                         Rotor.blade_data['c'][0], Rotor.blade_data['t/c'][0])
        self.assertEqual(tan_i, -16.635578670240793, 'should be -16.635578670240793 ')  # -16.635578670240793
        self.assertEqual(norm_i, 56.735026640634054, 'should be 56.735026640634054 ')  # 56.735026640634054

if __name__ == '__main__':
    unittest.main()
    rot = Rotor()
    Cl, Cd = rot.lift_drag_coeff(alpha=0, t_c=24.1)
    print(f'Cl {Cl}')
    print(f'Cd {Cd}')

    assert rot.lift_drag_coeff(alpha=0, t_c=24.1)[0] == 0.3391, 'review lift_drag_coeff method'
