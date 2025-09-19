import numpy as np


def estimate(idx, time, Dobs, t, P, V, B, f0, c=3e8):
    """
    Compute modeled Doppler (Dcom) at selected sample indices and return a comparison table.

    This is a python translation of the estimate.m function for argos simulator
    written by E. Schrama.

    Parameters
    ----------
    idx : (M,) array_like of int
        0-based indices selecting samples from the ephemeris arrays.
    time : (M,) array_like
        Observation times (s) corresponding to Dobs.
    Dobs : (M,) array_like
        Observed Doppler shifts (Hz).
    t : (N,) array_like
        Ephemeris epochs (s).
    P : (N,3) array_like
        Satellite inertial positions (m) at each epoch in t.
    V : (N,3) array_like
        Satellite inertial velocities (m/s) at each epoch in t.
    B : array_like
        Beacon inertial position(s) (m). Accepted shapes:
          - (3,)      : constant beacon position
          - (N,3)     : beacon position per ephemeris epoch (indexed by idx)
          - (M,3)     : beacon position per observation
    f0 : float
        Beacon transmit frequency (Hz).
    c : float, optional
        Speed of light (m/s). Default 3e8.

    Returns
    -------
    summary : (M, 5) ndarray
        Columns: [time, t[idx], Dobs, Dcom, (Dobs - Dcom)].
    Dcom : (M,) ndarray
        Modeled Doppler (Hz).

    Notes
    -----
    - `idx` is 0-based (Python convention).
    - All vectors are inertial and must be in consistent units/frames.
    """
    # --- Coerce & basic checks ---
    idx = np.asarray(idx, dtype=int).ravel()
    time = np.asarray(time, dtype=float).ravel()
    Dobs = np.asarray(Dobs, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    P = np.asarray(P, dtype=float)
    V = np.asarray(V, dtype=float)
    B = np.asarray(B, dtype=float)

    M = idx.size
    if not (time.size == M and Dobs.size == M):
        raise ValueError("idx, time, and Dobs must have the same length (M).")

    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("P must have shape (N,3).")
    if V.ndim != 2 or V.shape[1] != 3:
        raise ValueError("V must have shape (N,3).")
    N = P.shape[0]
    if t.size != N or V.shape[0] != N:
        raise ValueError("t, P, and V must share the same length N.")

    if np.any((idx < 0) | (idx >= N)):
        raise IndexError(
            "Some indices in 'idx' are out of bounds for the ephemeris arrays."
        )

    # --- Select satellite states at the observation indices ---
    Psel = P[idx]  # (M,3)
    Vsel = V[idx]  # (M,3)
    tsel = t[idx]  # (M,)

    # --- Prepare beacon positions aligned with observations ---
    if B.ndim == 1:
        if B.shape[0] != 3:
            raise ValueError("B with ndim==1 must have shape (3,).")
        Bsel = np.broadcast_to(B, (M, 3))
    elif B.ndim == 2 and B.shape[1] == 3:
        if B.shape[0] == N:
            Bsel = B[idx]  # per-epoch beacon, index by idx
        elif B.shape[0] == M:
            Bsel = B  # per-observation beacon
        else:
            raise ValueError("B must be (3,), (N,3), or (M,3).")
    else:
        raise ValueError("B must be (3,), (N,3), or (M,3).")

    # --- LOS, radial speed, Doppler ---
    # Line-of-sight vector from satellite to beacon: sat->beacon = P - B
    los = Psel - Bsel  # (M,3)
    dis = np.linalg.norm(los, axis=1)  # (M,)
    if np.any(dis == 0.0):
        raise ValueError(
            "Zero distance between satellite and beacon for some observation."
        )

    # Unit vector from satellite toward beacon: -los / |los|
    n_hat = -los / dis[:, None]  # (M,3)

    # Radial speed along LOS
    vee = np.sum(n_hat * Vsel, axis=1)  # (M,)

    # Modeled Doppler
    Dcom = (c / (c - vee)) * f0 - f0  # (M,)

    # Summary table
    summary = np.column_stack((time, tsel, Dobs, Dcom, Dobs - Dcom))
    return summary, Dcom
