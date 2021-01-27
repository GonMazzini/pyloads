import numpy as np


class BladeFeatures:
    # TODO:
    #   add functionality to check if specific blade data is correct
    def __init__(self, radio=np.array([2.8, 11.0027383, 16.871021, 22.9641225, 32.3076383,
                                       41.5730197, 50.4110617, 58.5344352, 65.7500952, 71.9674921,
                                       77.1859743, 78.7133469, 80.1402166, 82.7084849, 84.9251565,
                                       86.8264859, 88.4486629, 89.166])
                 , twist=np.array([14.5, 14.4320088, 12.5459743, 8.89245229, 6.38219779,
                                   4.67357583, 2.88655238, 1.20718981, -0.13037962, -1.11392478,
                                   -1.86195876, -2.0776268, -2.2800794, -2.64408409, -2.94507068,
                                   -3.17732896, -3.35752648, -3.428])
                 , cord=np.array([5.38, 5.45248001, 5.86600432, 6.18451124, 6.01996243,
                                  5.42189119, 4.69821952, 3.99867965, 3.39809449, 2.91304006,
                                  2.53594729, 2.4307714, 2.33143176, 2.12906893, 1.9033225,
                                  1.62601186, 1.18359898, 0.6])
                 , t_c=np.array([100., 86.0491476, 61.0969583, 43.0356446, 32.4150527,
                                 27.8103477, 25.317976, 24.2644205, 24.1006815, 24.1,
                                 24.1000016, 24.1000006, 24.1, 24.1, 24.1,
                                 24.1, 24.1, 24.1])):
        self.radio = radio
        self.twist = twist
        self.cord = cord
        self.t_c = t_c

    def plot_cord(self):
        import matplotlib.pyplot as plt
        plt.plot(self.radio, self.cord)
        plt.xlabel('radial position [m]')
        plt.ylabel('cord length [m]')
        plt.show()


if __name__ == '__main__':
    bld = BladeFeatures()
    print(bld.cord)
    bld.plot_cord()
