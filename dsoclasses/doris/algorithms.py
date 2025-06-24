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
