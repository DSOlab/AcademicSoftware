import numpy as np


def detect(idx, time, doppler, maximum=20.0):
    """
    Split observations into contiguous passes based on time gaps.

    This is a python translation of the detect.m function for argos simulator
    written by E. Schrama.

    Parameters
    ----------
    idx : array_like
        Unused; kept for API parity with MATLAB function.
    time : (N,) array_like
        Observation times (seconds). Assumed in chronological order.
    doppler : array_like
        Unused; kept for API parity with MATLAB function.
    maximum : float, optional
        Gap threshold (seconds). If time[i] - time[i-1] > maximum,
        a new pass starts. Default is 20.0.

    Returns
    -------
    npass : int
        Number of detected passes.
    pass_ids : (N,) ndarray of int
        Pass label for each observation, starting at 1.

    Notes
    -----
    - Behavior matches the MATLAB code: it does not sort by time.
      If `time` isn’t monotonically increasing, results may be odd.
    """
    t = np.asarray(time, dtype=float).ravel()
    n = t.size
    if n == 0:
        return 0, np.array([], dtype=int)

    pass_ids = np.empty(n, dtype=int)
    npass = 1
    pass_ids[0] = npass

    for i in range(1, n):
        if (t[i] - t[i - 1]) > maximum:
            npass += 1
        pass_ids[i] = npass

    return npass, pass_ids
