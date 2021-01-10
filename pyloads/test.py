from static_loads import Rotor

if __name__ == '__main__':
    rot=Rotor()
    Cl, Cd = rot.lift_drag_coeff(alpha=0, t_c=24.1)
    print(f'Cl {Cl}')
    print(f'Cd {Cd}')

    assert rot.lift_drag_coeff(alpha=0, t_c=24.1)[0] == 0.3391, 'review lift_drag_coeff method'




