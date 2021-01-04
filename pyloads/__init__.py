from static_loads import Rotor
import numpy as np

rot = Rotor()

print('given t/c(100%) and alpha (in degrees) calculate Cd and Cl')
print(rot.lift_drag_coeff(t=24.1, alpha=2*np.pi/180))
