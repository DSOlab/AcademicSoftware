import numpy as np

""" F0, aka USO frequency in [Hz]. """
USO_F0 = 5e6;

def beacon_nominal_frequency(shift_factor: int) -> tuple[float, float]:
    """ 
    Compute the S1 and U2 (aka 2 GHz and 400 MHz) nominal frequencies
    for a DORIS beacon.

    Parameters
    ----------
    shift_factor : int
        The beacon's shift factor (e.g., as extracted from the 'STATION REFERENCE'
        field from a DORIS RINEX file).

    Returns
    -------
    s1_freq : float
        The S1 (aka 2 GHz) nominal frequency [Hz].
    u2_freq : float
        The U2 (aka 400 MHz) nominal frequency [Hz].
    """
    two26 = 2**26
    fac1 = USO_F0 * 0.75e0
    fac2 = (USO_F0 * (87e0 * shift_factor)) / (5e0 * two26)
    s1_freq = 543e0 * fac1 + 543e0 * fac2
    u2_freq = 107e0 * fac1 + 107e0 * fac2
    return s1_freq, u2_freq

def starec_pcv(el):
    zen = 90. - np.degrees(el)
    angles = np.arange(0.0, 91.0, 5.0)  # 0,5,...,90  → length 19
    vals = np.array([0.00,    0.06,   -0.32,   -1.12,   -2.87,   -4.02,   -3.44,   -2.15,   -1.73,   -1.73,   -0.08,    1.37,    2.20,    5.37,    7.02,   10.70,  13.86,   17.27,   22.37])

    if vals.shape[0] != angles.size:
        raise ValueError(f"`vals` must have length {angles.size} (got {vals.shape[0]}).")

    a = np.asarray(zen, dtype=float)
    if np.any((a < 0) | (a > 90)):
        raise ValueError(f"Angle(s) `a` must be within [0, 90] degrees; give angle was {a:.2f} deg.")

    return np.interp(a, angles, vals)
