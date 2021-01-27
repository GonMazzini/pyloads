import pandas as pd


class Operation:

    def __init__(self):
        self.operation = pd.DataFrame(
            {'u': [4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.,
                   22., 23., 24., 25., ],
             'pitch': [2.751, 1.966, 0.896, 0., 0., 0., 0., 0., 4.502, 7.266,
                       9.292, 10.958, 12.499, 13.896, 15.2, 16.432, 17.618, 18.758, 19.86,
                       20.927,
                       21.963, 22.975],
             'RPM': [6., 6., 6., 6., 6.426, 7.229, 8.032, 8.836, 9.6, 9.6, 9.6, 9.6,
                     9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6, ]
             })

    def show_operation(self, u=30):
        """Show the pitch angle (degrees) and the RPM for different wind speeds.

         Parameters
        ----------
        u : int
        wind speed (should be between 4 and 25)
        """
        # TODO
        #   exception handling for u parameter
        if u == 30:

            out = self.operation  # DataFrame
        else:
            out = self.operation.iloc[u-4] # Tuple

        print(out)
        return out

if __name__ == "__main__":
    # WT_oper = pd.read_csv('operation.txt', sep='\s+')
    # WT_oper.index = WT_oper.u
    # print(np.array(WT_oper.RPM.values))
    op = Operation()
    u, pitch, rpm = op.show_operation(u=4)
