
---

## Satellite position at signal emission time

GNSS observations are time-tagged at the **receiver**: the RINEX epoch time is the *reception* time (t_r).
However, the signal left the satellite at an **earlier** time (t_e), because radio waves travel at finite speed (c).

Let

* (\mathbf{r}_r) = receiver position in ECEF (your `rsta`),
* (\mathbf{r}_s(t)) = satellite position in ECEF at time (t) (from SP3 interpolation),
* (c) = speed of light.

The basic **light-time equation** is

[
\big|\mathbf{r}_s(t_e) - \mathbf{r}_r\big|
= c,(t_r - t_e) \equiv c,\Delta t,
]

where

* (\Delta t = t_r - t_e) is the signal **time of flight**.

The reception time (t_r) is known; the emission time (t_e) is unknown and appears inside (\mathbf{r}_s(t_e)). Your `sat_at_emission_time` function solves this equation numerically.

---

## Algorithm: iterative light-time solution

We want to find (t_e) such that

[
\rho(t_e) = |\mathbf{r}_s(t_e) - \mathbf{r}_r| = c,(t_r - t_e).
]

Your code uses a simple fixed-point iteration:

1. **Initialization**

   Set the first guess of emission time equal to the reception time:
   [
   t_e^{(0)} = t_r.
   ]

2. **Iteration**

   For iteration (k = 0,1,2,\dots):

   1. Interpolate the satellite position at the current emission time guess:
      [
      \mathbf{r}_s^{(k)} = \mathbf{r}*s\big(t_e^{(k)}\big)
      = (x^{(k)}, y^{(k)}, z^{(k)}),
      ]
      and get the corresponding satellite clock correction (c*{\text{sat}}^{(k)}) from the SP3 interpolator.

   2. Compute the **geometric range** to the receiver (your `geometric_range`):
      [
      \rho^{(k)} = \big|\mathbf{r}_s^{(k)} - \mathbf{r}_r\big|.
      ]

   3. Compute the **signal travel time**:
      [
      \Delta t^{(k)} = \frac{\rho^{(k)}}{c}.
      ]

   4. Update the emission time:
      [
      t_e^{(k+1)} = t_r - \Delta t^{(k)}.
      ]

   5. Interpolate again at the new emission time:
      [
      \mathbf{r}_s^{(k+1)} = \mathbf{r}_s\big(t_e^{(k+1)}\big)
      = (x^{(k+1)}, y^{(k+1)}, z^{(k+1)}).
      ]

   6. Check convergence: in the code you compare the coordinate change between iterations against a **tolerance** vector (\boldsymbol{\tau}) (e.g. 0.1 m in each component). Formally:

      [
      \big|\mathbf{r}_s^{(k+1)} - \mathbf{r}_s^{(k)}\big|
      ===================================================

      \big(|x^{(k+1)} - x^{(k)}|,
      |y^{(k+1)} - y^{(k)}|,
      |z^{(k+1)} - z^{(k)}|\big)
      ]

      If **all components** are below the tolerance,
      [
      \big|\mathbf{r}_s^{(k+1)} - \mathbf{r}_s^{(k)}\big|
      < \boldsymbol{\tau},
      ]
      the iteration is considered converged.

      Otherwise, set (k \leftarrow k+1) and repeat, up to a maximum number of iterations.

3. **Outputs**

   Once convergence is reached, the function returns:

   * (x, y, z): the satellite ECEF coordinates at emission time (t_e);
   * (c_{\text{sat}}): the satellite clock correction at (t_e);
   * (\Delta t): the final signal travel time
     [
     \Delta t = \frac{\rho}{c} = t_r - t_e.
     ]

   So the function signature conceptually is:

   ```text
   (x_e, y_e, z_e, c_sat_e, Δt) = sat_at_emission_time(r_sta, t_r, ephemeris, sat_id)
   ```

   where (t_e = t_r - \Delta t).

---

## Why this correction is needed

If you naively used the satellite position at **reception time** (t_r) instead of emission time (t_e),

[
\mathbf{r}_s(t_r) \quad \text{instead of} \quad \mathbf{r}_s(t_e),
]

you would introduce a geometric error of order

[
\delta \rho \approx |\mathbf{v}_s| ,\Delta t,
]

where (\mathbf{v}_s) is the satellite velocity (about 3–4 km/s) and (\Delta t \sim 0.07)–0.08 s. This is typically on the order of **hundreds of meters** in range.

By iteratively solving for the emission time, your `sat_at_emission_time` function ensures that:

* the satellite coordinates used in the observation model are consistent with the actual signal path, and
* the geometric range (\rho) and the travel time (\Delta t) satisfy the light-time equation to within the specified tolerance.

This is essential for high-precision GNSS processing (PPP, geodesy, etc.), and is the reason you use `sat_at_emission_time` instead of a simple `sat_at(sat, t_reception)`.
