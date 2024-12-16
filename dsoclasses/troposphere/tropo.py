import numpy as np
from dsoclasses.geodesy import transformations
from dsoclasses.troposphere import gmf, gpt3

def tropo_delay(rsta, t, el, gpt3_grid):
    if not hasattr(tropo_delay, "car"):
        lat, lon, hgt = transformations.car2ell(rsta[0], rsta[1], rsta[2])
        tropo_delay.ell = np.array((lat, lon, hgt))
        tropo_delay.car = np.array(rsta)

    if np.any(np.greater_equal(np.abs(tropo_delay.car-np.array(rsta)), np.array((.1, .1, .1)))):
        lat, lon, hgt = transformations.car2ell(rsta[0], rsta[1], rsta[2])
        tropo_delay.ell = np.array((lat, lon, hgt))
        tropo_delay.car = np.array(rsta)

    lon, lat, hgt = tropo_delay.ell
    meteo = gpt3.gpt3(t, lon, lat, hgt, gpt3_grid)
    zhd = gpt3.saastamoinen_zhd(lat, hgt, meteo['p'])
    zwd = gpt3.askne_zwd(meteo['e'], meteo['Tm'], meteo['la'])
    gmfh, gmfw = gmf.gmf(t, lat, lon, hgt, np.pi/2-el)
    return zhd * gmfh + zwd * gmfw
