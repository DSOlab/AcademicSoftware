import numpy as np


def sph2xyz(ae, lat, lon):
    #
    # compute beacon coordinates in earth fixed frame
    #
    xb = ae * np.cos(lon) * np.cos(lat)
    yb = ae * np.sin(lon) * np.cos(lat)
    zb = ae * np.sin(lat)
    return np.array([xb, yb, zb])
