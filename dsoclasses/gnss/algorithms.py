import inspect
import numpy as np
from scipy.spatial.transform import Rotation as nRot
from dsoclasses.gnss import systems as gs
from dsoclasses.time.pyattotime import at2pt, fsec2asec


def geometric_range(rsta, rsat):
    """Return the geometric range, i.e. distance between two points.
    rsta and rsat can be lists or np.array
    """
    return np.linalg.norm(np.array(rsat) - np.array(rsta))


def sat_at_emission_time(
    rsta, t_reception, sp3_interpolator, sat, tolerance=(0.1, 0.1, 0.1)
):
    """Compute and return the satellite position at signal emission time.
    Parameters:
    rsta : Coordinates of observing station (ECEF, cartesian) in [m]
    t_reception: Time of reception, e.g. as recorded in RINEX
    interpolator: An Sp3Interpolator instance, where we can query sat's
           coordinates. I.e. `interpolator.sat_at(sat, t_emission)` should
           work.
    sat: The satellite's id, as in the interpolator (e.g. 'G03'). if the
    interpolator only has one satellite, this is irrelevant.
    tolerance: An array, or list of max discrepancies untill we consider
         successeful convergence. These are actually maximum allowed
         coordinate discrepancies (in [m]).
    Returns:
    x, y, z, c, dt
        where x, y, z, are the (ECEF, cartesian) coordinates of the
    satellite at emission time; c is the satellite's clock correction at
    the same epoch. dt is t_reception - t_emission in fractional seconds.
    """
    MAX_ITERATIONS = 10
    # starting t_(emission) is t_(reception)
    t_emission = t_reception
    # get satellite coordinates at t, using the interpolator
    # x, y, z, clk = interpolator.sat_at(sat, t_emission)
    x, y, z, c = sp3_interpolator.sat_at(sat, t_emission)

    count_it = 0
    while count_it < MAX_ITERATIONS:
        # iterate while the is no 'significant' change in sat coordinates
        r = geometric_range((x, y, z), rsta)
        dt = r / gs.C
        t_emission = t_reception - fsec2asec(dt)
        xnew, ynew, znew, clk = sp3_interpolator.sat_at(sat, t_emission)
        if np.any(
            np.greater_equal(
                np.abs(np.array((xnew, ynew, znew)) - np.array((x, y, z))),
                np.array(tolerance),
            )
        ):
            x = xnew
            y = ynew
            z = znew
            count_it += 1
        else:
            return xnew, ynew, znew, clk, dt

    raise RuntimeError("[ERROR] Failed finding emission time for sat {:}".format(sat))
