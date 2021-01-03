from static_loads import Rotor
import numpy as np

rot = Rotor(rho=1.225)
a=2.8*np.pi/180
print(rot.lift_drag_coeff(t=5.38000000E+00, alpha=a))
