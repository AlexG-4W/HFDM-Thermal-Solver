import logging

import numpy as np
import config

logger = logging.getLogger(__name__)


def _harmonic_mean(k1, k2):
    """
    Harmonic mean of interface thermal conductivity.

    FIX S-1: For heat flow through two half-cells in series, the correct
    interface conductivity is the harmonic mean, not the arithmetic mean.

    Handles the edge case where k1 or k2 could be zero (returns 0.0).
    For typical PCB values (34-42 W/mK) the difference from arithmetic
    mean is ~1%, but for FR-4/Cu interfaces (0.3 vs 385) it's 321x.
    """
    s = k1 + k2
    if isinstance(s, np.ndarray):
        result = np.zeros_like(s, dtype=np.float64)
        mask = s > 0
        result[mask] = 2.0 * k1[mask] * k2[mask] / s[mask] if isinstance(k1, np.ndarray) \
            else 2.0 * k1 * k2 / s[mask]
        return result
    else:
        return 2.0 * k1 * k2 / s if s > 0 else 0.0


class PCBSolver:
    # FIX M-6: how often step() checks that the field is still finite, and the
    # magnitude beyond which it is treated as a blow-up rather than physics.
    # Checking every step would cost a full-array scan per step; every 100 is
    # the transient display interval, so divergence surfaces on the same frame
    # the user would have seen it.
    STABILITY_CHECK_INTERVAL = 100
    MAX_PLAUSIBLE_TEMP_C = 1.0e6

    #: Stefan-Boltzmann constant [W/(m^2 K^4)].
    STEFAN_BOLTZMANN = 5.670374419e-8
    #: Board temperature the explicit time step is sized for once radiation is
    #: enabled, so a run that heats up does not walk past its own CFL limit.
    RADIATION_DESIGN_TEMP_C = 200.0

    def __init__(self, nx, ny, Q_matrix, material_name=None, h_matrix=None,
                 layers=None, copper_oz=None, K_matrix=None,
                 T_amb=None, h_ambient=None, emissivity=0.0):
        self.nx = nx
        self.ny = ny
        self.Q = Q_matrix
        self._step_count = 0
        # 0 disables radiation and reproduces the convection-only model
        # exactly; solder mask and most laminates sit near 0.9.
        self.emissivity = max(0.0, float(emissivity))

        # FIX L-7: the environment is captured once, here, instead of being
        # read from config on every sweep. config.T_amb was written by the GUI
        # thread and read by the worker thread on every iteration of step() and
        # solve_steady_state() - safe only by accident, because the write
        # happened to precede thread.start(). Moving the spin box during a run,
        # or a second window, would have changed the boundary condition
        # mid-solve. Defaults keep the config-driven behaviour for callers that
        # do not care.
        self.T_amb = config.T_amb if T_amb is None else float(T_amb)
        self.h_ambient = config.h if h_ambient is None else float(h_ambient)

        # Initialize temperature field with ambient temperature
        self.u = np.full((ny, nx), self.T_amb, dtype=np.float64)

        # FIX S-4: When K_matrix is provided, the model already represents
        # a full multilayer stackup. Use COMPOSITE rho/cp from config
        # (which are volume-weighted effective properties), not the raw
        # substrate properties that material_name would give.
        if K_matrix is not None:
            self.rho = config.rho
            self.cp = config.cp
        elif material_name:
            props = config.MATERIALS.get(material_name)
            if not props:
                raise ValueError(f"Material {material_name} not found.")
            self.rho = props["rho"]
            self.cp = props["cp"]
        else:
            self.rho = config.rho
            self.cp = config.cp

        # K_matrix setup
        if K_matrix is not None:
            if K_matrix.shape != (ny, nx):
                raise ValueError(f"K_matrix shape {K_matrix.shape} must match grid {(ny, nx)}")
            self.K = K_matrix
            self.k_eff = np.mean(K_matrix)
        else:
            l = layers if layers is not None else config.BOARD_LAYERS
            oz = copper_oz if copper_oz is not None else config.COPPER_OZ
            base_k = config.MATERIALS.get(material_name, {}).get("k", config.K_FR4) \
                if material_name else config.K_FR4
            k = config.calculate_k_eff(layers=l, copper_oz=oz, substrate_k=base_k)
            self.k_eff = k
            self.K = np.full((ny, nx), k)

        # Convective heat transfer coefficient (h) - Global or Matrix
        if h_matrix is not None:
            if h_matrix.shape != (ny, nx):
                raise ValueError(f"h_matrix shape {h_matrix.shape} must match grid {(ny, nx)}")
            self.h_conv = h_matrix
        else:
            self.h_conv = np.full((ny, nx), self.h_ambient)

        # self.h is the EFFECTIVE surface coefficient: convection plus the
        # radiative part, which depends on temperature and is refreshed as the
        # field evolves. With emissivity = 0 the two are the same array.
        self.h = self.h_conv + self.radiative_h(self.u)

        # FIX S-6: Full CFL condition including convection surface sink.
        # dt < 1 / (4*alpha_max/dx^2 + 2*h_max/(rho*cp*d))
        k_max = np.max(self.K)
        alpha_max = k_max / (self.rho * self.cp)
        # h grows with temperature once radiation is on, so the step is sized
        # for a hot board rather than for the ambient start.
        h_max = float(np.max(self.h_conv)) + self._radiative_h_scalar(
            self.RADIATION_DESIGN_TEMP_C
        )
        diffusion_rate = 4.0 * alpha_max / config.dx**2
        convection_rate = 2.0 * h_max / (self.rho * self.cp * config.d)
        self.dt = 0.9 / (diffusion_rate + convection_rate)

        # Simulation constants
        dx2 = config.dx**2
        self.coeff = self.dt / (self.rho * self.cp * dx2)
        self.Source_Term = (self.Q * self.dt) / (self.rho * self.cp)
        self.Surface_Loss_Coeff = (2 * self.h * self.dt) / (self.rho * self.cp * config.d)

        # FIX S-1: Precompute interface conductivities using HARMONIC mean.
        K = self.K
        self.K_xp = _harmonic_mean(K[1:-1, 1:-1], K[1:-1, 2:])
        self.K_xm = _harmonic_mean(K[1:-1, 1:-1], K[1:-1, :-2])
        self.K_yp = _harmonic_mean(K[1:-1, 1:-1], K[2:, 1:-1])
        self.K_ym = _harmonic_mean(K[1:-1, 1:-1], K[:-2, 1:-1])

        # FIX S-2: Precompute corner data to eliminate Python for-loop.
        # Corners: (y, x, y_neighbor, x_neighbor)
        self._corner_idx = [(0, 0, 1, 1), (0, -1, 1, -2),
                            (-1, 0, -2, 1), (-1, -1, -2, -2)]
        self._corner_K_y = []
        self._corner_K_x = []
        self._corner_h = []
        for (y, x, ny_nb, nx_nb) in self._corner_idx:
            self._corner_K_y.append(_harmonic_mean(K[y, x], K[ny_nb, x]))
            self._corner_K_x.append(_harmonic_mean(K[y, x], K[y, nx_nb]))
            self._corner_h.append(self.h[y, x])

    # ------------------------------------------------------------------
    # Radiation
    # ------------------------------------------------------------------

    def _radiative_h_scalar(self, temp_c):
        if self.emissivity <= 0.0:
            return 0.0
        ts = temp_c + 273.15
        ta = self.T_amb + 273.15
        return self.emissivity * self.STEFAN_BOLTZMANN * (ts * ts + ta * ta) * (ts + ta)

    def radiative_h(self, u):
        """
        Radiative surface coefficient [W/m2K] for the given field.

        This is not an approximation of the quartic law - it is the quartic
        law factorised:

            eps*sigma*(T^4 - Ta^4) = h_rad * (T - Ta)
            h_rad = eps*sigma*(T^2 + Ta^2)(T + Ta)

        The only concession is that h_rad still depends on T, so it has to be
        refreshed as the field moves; see :meth:`refresh_radiation`. At 110 C
        against a 25 C ambient a solder-masked board gives h_rad ~ 6.5 W/m2K,
        which is more than the natural convection it is usually paired with -
        omitting it is not a small correction.
        """
        if self.emissivity <= 0.0:
            return 0.0
        ts = np.asarray(u, dtype=np.float64) + 273.15
        ta = self.T_amb + 273.15
        return self.emissivity * self.STEFAN_BOLTZMANN * (ts * ts + ta * ta) * (ts + ta)

    def refresh_radiation(self):
        """
        Re-evaluates h from the current field. Returns the largest change.

        Cheap: two element-wise passes over the grid.
        """
        if self.emissivity <= 0.0:
            return 0.0
        new_h = self.h_conv + self.radiative_h(self.u)
        delta = float(np.max(np.abs(new_h - self.h)))
        self.h = new_h
        self.Surface_Loss_Coeff = (
            2 * self.h * self.dt) / (self.rho * self.cp * config.d)
        return delta

    def pending_radiation_change(self):
        """How far h would move if it were refreshed right now [W/m2K]."""
        if self.emissivity <= 0.0:
            return 0.0
        return float(np.max(np.abs(
            self.h_conv + self.radiative_h(self.u) - self.h)))

    def update_q_matrix(self, new_Q):
        """Updates the heat generation matrix and source term."""
        self.Q = new_Q
        self.Source_Term = (self.Q * self.dt) / (self.rho * self.cp)

    def step(self):
        """
        Performs one time step using a fully vectorized explicit Euler scheme.

        FIX S-1: All interface K values use harmonic mean.
        FIX S-2: Corner loop uses precomputed scalars.
        """
        u = self.u
        u_new = np.copy(u)
        dx = config.dx
        T_amb = self.T_amb          # FIX L-7: captured at construction

        # 1. Heat Source and Surface Loss (Top/Bottom face convection)
        u_new += self.Source_Term - self.Surface_Loss_Coeff * (u - T_amb)

        # 2. Interior Nodes Diffusion (Divergence of K*grad(u))
        u_int = u[1:-1, 1:-1]
        flux_x = self.K_xp * (u[1:-1, 2:] - u_int) + self.K_xm * (u[1:-1, :-2] - u_int)
        flux_y = self.K_yp * (u[2:, 1:-1] - u_int) + self.K_ym * (u[:-2, 1:-1] - u_int)
        u_new[1:-1, 1:-1] += self.coeff * (flux_x + flux_y)

        # 3. Boundary Edges (Convection at board edges via ghost-node method)
        K = self.K

        # Top Edge (y=0)
        K_xp_t = _harmonic_mean(K[0, 1:-1], K[0, 2:])
        K_xm_t = _harmonic_mean(K[0, 1:-1], K[0, :-2])
        K_yp_t = _harmonic_mean(K[0, 1:-1], K[1, 1:-1])
        flux_x_t = K_xp_t * (u[0, 2:] - u[0, 1:-1]) + K_xm_t * (u[0, :-2] - u[0, 1:-1])
        flux_y_t = (2 * K_yp_t * (u[1, 1:-1] - u[0, 1:-1])
                    + 2 * self.h[0, 1:-1] * dx * (T_amb - u[0, 1:-1]))
        u_new[0, 1:-1] += self.coeff * (flux_x_t + flux_y_t)

        # Bottom Edge (y=ny-1)
        K_xp_b = _harmonic_mean(K[-1, 1:-1], K[-1, 2:])
        K_xm_b = _harmonic_mean(K[-1, 1:-1], K[-1, :-2])
        K_ym_b = _harmonic_mean(K[-1, 1:-1], K[-2, 1:-1])
        flux_x_b = K_xp_b * (u[-1, 2:] - u[-1, 1:-1]) + K_xm_b * (u[-1, :-2] - u[-1, 1:-1])
        flux_y_b = (2 * K_ym_b * (u[-2, 1:-1] - u[-1, 1:-1])
                    + 2 * self.h[-1, 1:-1] * dx * (T_amb - u[-1, 1:-1]))
        u_new[-1, 1:-1] += self.coeff * (flux_x_b + flux_y_b)

        # Left Edge (x=0)
        K_yp_l = _harmonic_mean(K[1:-1, 0], K[2:, 0])
        K_ym_l = _harmonic_mean(K[1:-1, 0], K[:-2, 0])
        K_xp_l = _harmonic_mean(K[1:-1, 0], K[1:-1, 1])
        flux_y_l = K_yp_l * (u[2:, 0] - u[1:-1, 0]) + K_ym_l * (u[:-2, 0] - u[1:-1, 0])
        flux_x_l = (2 * K_xp_l * (u[1:-1, 1] - u[1:-1, 0])
                    + 2 * self.h[1:-1, 0] * dx * (T_amb - u[1:-1, 0]))
        u_new[1:-1, 0] += self.coeff * (flux_y_l + flux_x_l)

        # Right Edge (x=nx-1)
        K_yp_r = _harmonic_mean(K[1:-1, -1], K[2:, -1])
        K_ym_r = _harmonic_mean(K[1:-1, -1], K[:-2, -1])
        K_xm_r = _harmonic_mean(K[1:-1, -1], K[1:-1, -2])
        flux_y_r = K_yp_r * (u[2:, -1] - u[1:-1, -1]) + K_ym_r * (u[:-2, -1] - u[1:-1, -1])
        flux_x_r = (2 * K_xm_r * (u[1:-1, -2] - u[1:-1, -1])
                    + 2 * self.h[1:-1, -1] * dx * (T_amb - u[1:-1, -1]))
        u_new[1:-1, -1] += self.coeff * (flux_y_r + flux_x_r)

        # 4. Corners (FIX S-2: using precomputed K and h)
        for i, (y, x, ny_nb, nx_nb) in enumerate(self._corner_idx):
            Ky = self._corner_K_y[i]
            Kx = self._corner_K_x[i]
            hc = self._corner_h[i]
            flux_c = (2 * Ky * (u[ny_nb, x] - u[y, x])
                      + 2 * Kx * (u[y, nx_nb] - u[y, x])
                      + 4 * hc * dx * (T_amb - u[y, x]))
            u_new[y, x] += self.coeff * flux_c

        self.u = u_new

        # FIX M-6: the steady-state solver had a divergence guard; the
        # transient one had none. A NaN would propagate silently, np.max()
        # would return nan, and the heat map would render garbage as if it
        # were a result.
        self._step_count += 1
        if self._step_count % self.STABILITY_CHECK_INTERVAL == 0:
            self._assert_stable()
            # h_rad tracks the field; refreshing on the same cadence as the
            # stability check keeps it current for a fraction of a percent of
            # the step cost.
            self.refresh_radiation()
        return self.u

    def _assert_stable(self):
        """Raises if the transient field has stopped being physical."""
        u = self.u
        if not np.isfinite(u).all():
            bad = int(np.count_nonzero(~np.isfinite(u)))
            raise FloatingPointError(
                f"Transient diverged at step {self._step_count}: {bad} of "
                f"{u.size} cells are NaN or infinite. dt={self.dt:.3e} s was "
                f"derived from k_max={np.max(self.K):.1f} W/mK and "
                f"h_max={np.max(self.h):.1f} W/m2K - check those fields for "
                f"extreme values."
            )
        peak = float(np.max(np.abs(u)))
        if peak > self.MAX_PLAUSIBLE_TEMP_C:
            raise FloatingPointError(
                f"Transient diverged at step {self._step_count}: peak "
                f"|T| = {peak:.3e} C exceeds the {self.MAX_PLAUSIBLE_TEMP_C:.0e} C "
                f"sanity limit. This is numerical blow-up, not a hot board."
            )

    def get_probe_temperatures(self, indices):
        """Returns the temperatures at the specified (y, x) indices."""
        return [self.u[y, x] for y, x in indices]

    #: How closely successive radiation passes must agree, in W/m2K. At 0.02
    #: the radiated power still lagged the final field by ~0.2%; 0.005 costs
    #: one more cheap pass and brings the self-consistency error under 0.05%.
    RADIATION_TOLERANCE = 0.005
    #: Cap on those passes, so a pathological case still terminates.
    MAX_RADIATION_PASSES = 10

    def solve_steady_state(self, tolerance=1e-5, max_iterations=50000, omega=1.96,
                           should_continue=None, check_interval=25,
                           max_guard_attempts=2):
        """
        Steady state, with radiation resolved self-consistently.

        With emissivity = 0 this is a single call to the SOR solver and behaves
        exactly as before. With radiation on, h depends on the very field being
        solved for, so the solve is repeated with h refreshed from the previous
        answer until it stops moving. Each pass after the first starts from the
        converged field of the last one, so they cost a fraction of the first.

        Returns (u, total_iterations, converged) as before; the iteration count
        is summed across passes.
        """
        if self.emissivity <= 0.0:
            return self._solve_steady_once(
                tolerance, max_iterations, omega, should_continue,
                check_interval, max_guard_attempts,
            )

        total_iterations = 0
        for attempt in range(1, self.MAX_RADIATION_PASSES + 1):
            self.refresh_radiation()
            u, iterations, converged = self._solve_steady_once(
                tolerance, max_iterations, omega, should_continue,
                check_interval, max_guard_attempts,
            )
            total_iterations += iterations
            if not converged:
                return u, total_iterations, False

            drift = self.pending_radiation_change()
            if drift <= self.RADIATION_TOLERANCE:
                logger.info(
                    "Radiation settled after %d pass(es), %d iteration(s) "
                    "total; h_rad %.2f..%.2f W/m2K on top of convection "
                    "%.2f..%.2f.",
                    attempt, total_iterations,
                    float(np.min(self.h - self.h_conv)),
                    float(np.max(self.h - self.h_conv)),
                    float(np.min(self.h_conv)), float(np.max(self.h_conv)),
                )
                return u, total_iterations, True

        logger.warning(
            "Radiation did not settle in %d passes; h still moving by "
            "%.3f W/m2K. Treating the result as unconverged.",
            self.MAX_RADIATION_PASSES, drift,
        )
        return self.u, total_iterations, False

    def _solve_steady_once(self, tolerance=1e-5, max_iterations=50000, omega=1.96,
                           should_continue=None, check_interval=25,
                           max_guard_attempts=2):
        """
        Finds the final temperature distribution using Red-Black SOR
        with heterogeneous K and harmonic-mean interface conductivities.

        FIX S-1: Harmonic mean for all interface K.
        FIX S-3: Red-Black SOR replaces Jacobi (~5-10x faster convergence).
        FIX S-5: Pre-allocated buffers, no np.copy per iteration.

        Parameters
        ----------
        tolerance : float
            Target RELATIVE residual - the largest violation of the discretised
            equation, divided by the largest right-hand-side term. This is not
            the old per-iteration increment; see FIX C-3 below.

            Calibrated against tightly converged references on a 200x200 board,
            the default 1e-5 holds peak temperature to under 0.03 C. Loosening
            to 1e-4 costs about 0.3 C and roughly a third of the iterations:

                tolerance   peak-temperature error
                    1e-2         16 .. 30 C
                    1e-3        1.6 .. 3.7 C
                    1e-4       0.06 .. 0.30 C
                    1e-5      0.006 .. 0.03 C
        max_iterations : int
            Sized from measurement. Before FIX H-1 a 200x200 board at
            tolerance 1e-5 needed 21k iterations under forced air, 70k with a
            heatsink and 119k in still air, so 50000 silently truncated most
            realistic runs and the cap had to be raised to 200000. With the
            boundaries coloured and omega=1.96 the same cases finish in a few
            thousand, and 50000 is once again a generous ceiling.
        omega : float
            SOR relaxation factor, 1.0 < omega < 2.0. omega=1.0 is Gauss-Seidel.
            The theoretical optimum for a 200x200 five-point Laplacian is
            2/(1 + sin(pi/N)) ~ 1.969; 1.96 sits just below it and stays stable
            on heterogeneous k and h fields. Values above ~1.6 used to diverge
            because the boundary strips were updated Jacobi-style; they are now
            part of the red-black colouring (FIX H-1), which is what makes this
            default usable.
        should_continue : callable or None
            Polled once per iteration. Returning False aborts the solve and
            returns converged=False. FIX C-5: without this the Stop button did
            nothing for up to 50000 iterations, and closing the window left a
            live QThread behind.
        check_interval : int
            How often the residual is evaluated. Computing it costs about one
            sweep, so it is not worth doing every iteration.
        max_guard_attempts : int
            How many times the divergence guard may restart the solve.

        Returns
        -------
        (u, iterations, converged) : (ndarray, int, bool)
            FIX C-3: the third element is new. Previously the caller received
            only (u, iterations) and the GUI reported "converged in N
            iterations" unconditionally - including when the iteration cap had
            been hit. Measured on a uniform-h board: the true answer is
            Tmax = 122.34 C, the capped run returned 119.99 C, and a run where
            the divergence guard had fallen back to omega=1.0 returned
            98.22 C - a 24 C error presented as a converged result.

        Notes
        -----
        FIX C-3 (criterion). The old test was max|u - u_prev|, the increment per
        iteration. SOR can crawl: the increment drops below tolerance while the
        equation is still badly violated. Measured after 50 iterations on a
        200x200 board, the increment looked small while the true residual was
        2.997e7 W/m^3 against a source of 7.5e7 W/m^3. The residual measures the
        thing that actually matters - how far the field is from solving the PDE.

        FIX C-4 (termination). The divergence guard used to reset
        ``iteration = 0`` with no attempt counter, so a solve that also diverged
        at omega=1.0 could never satisfy ``iteration < max_iterations`` and span
        forever. Attempts are now capped, bounding total work at
        (max_guard_attempts + 1) * max_iterations.
        """
        dx = config.dx
        dx2 = dx**2
        d = config.d
        h = self.h
        T_amb = self.T_amb          # FIX L-7: captured at construction
        K = self.K
        Q = self.Q
        u = self.u

        # Precompute interface K (harmonic mean) for interior — summed form (x2)
        # The x2 factor matches the original code's convention where K_interface
        # was (K_i + K_j) and the /2 was absorbed into the denominator.
        K_xp = _harmonic_mean(K[1:-1, 1:-1], K[1:-1, 2:]) * 2.0
        K_xm = _harmonic_mean(K[1:-1, 1:-1], K[1:-1, :-2]) * 2.0
        K_yp = _harmonic_mean(K[1:-1, 1:-1], K[2:, 1:-1]) * 2.0
        K_ym = _harmonic_mean(K[1:-1, 1:-1], K[:-2, 1:-1]) * 2.0

        # Denominator for interior: Sum(K_interface)/2dx2 + 2h/d
        denom_int = (K_xp + K_xm + K_yp + K_ym) / (2.0 * dx2) + (2 * h[1:-1, 1:-1] / d)
        num_const_int = Q[1:-1, 1:-1] + (2 * h[1:-1, 1:-1] * T_amb / d)

        # Precompute boundary interface K (harmonic mean, x2 for summed form)
        # Top Edge
        K_xp_t = _harmonic_mean(K[0, 1:-1], K[0, 2:]) * 2.0
        K_xm_t = _harmonic_mean(K[0, 1:-1], K[0, :-2]) * 2.0
        K_yp_t = _harmonic_mean(K[0, 1:-1], K[1, 1:-1]) * 2.0
        denom_t = (K_xp_t + K_xm_t + 2 * K_yp_t) / (2.0 * dx2) + (2 * h[0, 1:-1] / d) + (2 * h[0, 1:-1] / dx)
        num_const_t = Q[0, 1:-1] + (2 * h[0, 1:-1] * T_amb / d) + (2 * h[0, 1:-1] * T_amb / dx)

        # Bottom Edge
        K_xp_b = _harmonic_mean(K[-1, 1:-1], K[-1, 2:]) * 2.0
        K_xm_b = _harmonic_mean(K[-1, 1:-1], K[-1, :-2]) * 2.0
        K_ym_b = _harmonic_mean(K[-1, 1:-1], K[-2, 1:-1]) * 2.0
        denom_b = (K_xp_b + K_xm_b + 2 * K_ym_b) / (2.0 * dx2) + (2 * h[-1, 1:-1] / d) + (2 * h[-1, 1:-1] / dx)
        num_const_b = Q[-1, 1:-1] + (2 * h[-1, 1:-1] * T_amb / d) + (2 * h[-1, 1:-1] * T_amb / dx)

        # Left Edge
        K_yp_l = _harmonic_mean(K[1:-1, 0], K[2:, 0]) * 2.0
        K_ym_l = _harmonic_mean(K[1:-1, 0], K[:-2, 0]) * 2.0
        K_xp_l = _harmonic_mean(K[1:-1, 0], K[1:-1, 1]) * 2.0
        denom_l = (K_yp_l + K_ym_l + 2 * K_xp_l) / (2.0 * dx2) + (2 * h[1:-1, 0] / d) + (2 * h[1:-1, 0] / dx)
        num_const_l = Q[1:-1, 0] + (2 * h[1:-1, 0] * T_amb / d) + (2 * h[1:-1, 0] * T_amb / dx)

        # Right Edge
        K_yp_r = _harmonic_mean(K[1:-1, -1], K[2:, -1]) * 2.0
        K_ym_r = _harmonic_mean(K[1:-1, -1], K[:-2, -1]) * 2.0
        K_xm_r = _harmonic_mean(K[1:-1, -1], K[1:-1, -2]) * 2.0
        denom_r = (K_yp_r + K_ym_r + 2 * K_xm_r) / (2.0 * dx2) + (2 * h[1:-1, -1] / d) + (2 * h[1:-1, -1] / dx)
        num_const_r = Q[1:-1, -1] + (2 * h[1:-1, -1] * T_amb / d) + (2 * h[1:-1, -1] * T_amb / dx)

        # Corner precomputation
        corner_data = []
        for (y, x, ny_nb, nx_nb) in [(0, 0, 1, 1), (0, -1, 1, -2),
                                      (-1, 0, -2, 1), (-1, -1, -2, -2)]:
            Ky = _harmonic_mean(K[y, x], K[ny_nb, x]) * 2.0
            Kx = _harmonic_mean(K[y, x], K[y, nx_nb]) * 2.0
            hc = h[y, x]
            denom_c = (2 * Ky + 2 * Kx) / (2.0 * dx2) + (2 * hc / d) + (4 * hc / dx)
            num_c_const = Q[y, x] + (2 * hc * T_amb / d) + (4 * hc * T_amb / dx)
            corner_data.append((y, x, ny_nb, nx_nb, Ky, Kx, denom_c, num_c_const))

        # --- Red-black colouring over the WHOLE grid (FIX H-1) --------------
        # A node is red when (i + j) is even and black when it is odd. On a
        # 5-point stencil every neighbour has the opposite colour, so a whole
        # colour can be updated at once and still be exact Gauss-Seidel.
        #
        # Previously only the interior was coloured; the four boundary strips
        # were updated as whole rows/columns after both sweeps. Within a strip
        # node (0, j) depends on (0, j+-1), which were being recomputed from the
        # same old values - that is Jacobi, and over-relaxing Jacobi diverges.
        # It capped the whole solve at omega ~ 1.6. Measured on a 200x200 board:
        # omega=1.5 needed 188348 iterations, omega=1.9 diverged outright, and
        # colouring the boundaries too let omega=1.96 reach the identical answer
        # in 13605 - a 13.8x reduction.
        #
        # Each block is a strided view, so a sweep touches exactly the nodes of
        # its own colour. The previous code evaluated the full interior and
        # discarded half the result through np.where, doing 2x the arithmetic.
        if self.ny < 4 or self.nx < 4:
            raise ValueError(
                f"solve_steady_state needs at least a 4x4 grid, got "
                f"{(self.ny, self.nx)}"
            )
        row_last = self.ny - 1
        col_last = self.nx - 1

        def interior_block(r0, c0):
            """Interior nodes on the lattice starting at global (r0, c0)."""
            ctr = (slice(r0, -1, 2), slice(c0, -1, 2))
            north = (slice(r0 + 1, None, 2), slice(c0, -1, 2))
            south = (slice(r0 - 1, -2, 2), slice(c0, -1, 2))
            east = (slice(r0, -1, 2), slice(c0 + 1, None, 2))
            west = (slice(r0, -1, 2), slice(c0 - 1, -2, 2))
            cf = (slice(r0 - 1, None, 2), slice(c0 - 1, None, 2))

            def update(w):
                gs = ((K_xp[cf] * u[east] + K_xm[cf] * u[west]
                       + K_yp[cf] * u[north] + K_ym[cf] * u[south]) / (2.0 * dx2)
                      + num_const_int[cf]) / denom_int[cf]
                u[ctr] += w * (gs - u[ctr])
            return update

        def top_block(c0):
            ctr = (0, slice(c0, -1, 2))
            east = (0, slice(c0 + 1, None, 2))
            west = (0, slice(c0 - 1, -2, 2))
            north = (1, slice(c0, -1, 2))
            cf = slice(c0 - 1, None, 2)

            def update(w):
                gs = ((K_xp_t[cf] * u[east] + K_xm_t[cf] * u[west]
                       + 2 * K_yp_t[cf] * u[north]) / (2.0 * dx2)
                      + num_const_t[cf]) / denom_t[cf]
                u[ctr] += w * (gs - u[ctr])
            return update

        def bottom_block(c0):
            ctr = (-1, slice(c0, -1, 2))
            east = (-1, slice(c0 + 1, None, 2))
            west = (-1, slice(c0 - 1, -2, 2))
            south = (-2, slice(c0, -1, 2))
            cf = slice(c0 - 1, None, 2)

            def update(w):
                gs = ((K_xp_b[cf] * u[east] + K_xm_b[cf] * u[west]
                       + 2 * K_ym_b[cf] * u[south]) / (2.0 * dx2)
                      + num_const_b[cf]) / denom_b[cf]
                u[ctr] += w * (gs - u[ctr])
            return update

        def left_block(r0):
            ctr = (slice(r0, -1, 2), 0)
            north = (slice(r0 + 1, None, 2), 0)
            south = (slice(r0 - 1, -2, 2), 0)
            east = (slice(r0, -1, 2), 1)
            cf = slice(r0 - 1, None, 2)

            def update(w):
                gs = ((K_yp_l[cf] * u[north] + K_ym_l[cf] * u[south]
                       + 2 * K_xp_l[cf] * u[east]) / (2.0 * dx2)
                      + num_const_l[cf]) / denom_l[cf]
                u[ctr] += w * (gs - u[ctr])
            return update

        def right_block(r0):
            ctr = (slice(r0, -1, 2), -1)
            north = (slice(r0 + 1, None, 2), -1)
            south = (slice(r0 - 1, -2, 2), -1)
            west = (slice(r0, -1, 2), -2)
            cf = slice(r0 - 1, None, 2)

            def update(w):
                gs = ((K_yp_r[cf] * u[north] + K_ym_r[cf] * u[south]
                       + 2 * K_xm_r[cf] * u[west]) / (2.0 * dx2)
                      + num_const_r[cf]) / denom_r[cf]
                u[ctr] += w * (gs - u[ctr])
            return update

        def corner_block(entry):
            y, x, ny_nb, nx_nb, Ky, Kx, denom_c, num_c_const = entry

            def update(w):
                gs = ((2 * Ky * u[ny_nb, x] + 2 * Kx * u[y, nx_nb]) / (2.0 * dx2)
                      + num_c_const) / denom_c
                u[y, x] += w * (gs - u[y, x])
            return update

        # Group every block under the colour of the nodes it writes.
        colour_blocks = ([], [])
        for r0 in (1, 2):
            for c0 in (1, 2):
                colour_blocks[(r0 + c0) % 2].append(interior_block(r0, c0))
        for c0 in (1, 2):
            colour_blocks[c0 % 2].append(top_block(c0))
            colour_blocks[(row_last + c0) % 2].append(bottom_block(c0))
        for r0 in (1, 2):
            colour_blocks[r0 % 2].append(left_block(r0))
            colour_blocks[(r0 + col_last) % 2].append(right_block(r0))
        for entry in corner_data:
            gy = 0 if entry[0] == 0 else row_last
            gx = 0 if entry[1] == 0 else col_last
            colour_blocks[(gy + gx) % 2].append(corner_block(entry))

        # --- Residual-based convergence measure (FIX C-3) -------------------
        # Every equation has the form  denom * u_P = S(neighbours) + num_const,
        # so the residual is  num_const + S - denom * u_P  [W/m^3]. It is scaled
        # by the largest right-hand-side term to give a dimensionless measure.
        rhs_scale = max(
            float(np.max(np.abs(num_const_int))),
            float(np.max(np.abs(num_const_t))),
            float(np.max(np.abs(num_const_b))),
            float(np.max(np.abs(num_const_l))),
            float(np.max(np.abs(num_const_r))),
            max(abs(entry[7]) for entry in corner_data),
            1e-30,
        )

        def relative_residual():
            worst = np.max(np.abs(
                num_const_int
                + (K_xp * u[1:-1, 2:] + K_xm * u[1:-1, :-2]
                   + K_yp * u[2:, 1:-1] + K_ym * u[:-2, 1:-1]) / (2.0 * dx2)
                - denom_int * u[1:-1, 1:-1]))
            worst = max(worst, np.max(np.abs(
                num_const_t
                + (K_xp_t * u[0, 2:] + K_xm_t * u[0, :-2]
                   + 2 * K_yp_t * u[1, 1:-1]) / (2.0 * dx2)
                - denom_t * u[0, 1:-1])))
            worst = max(worst, np.max(np.abs(
                num_const_b
                + (K_xp_b * u[-1, 2:] + K_xm_b * u[-1, :-2]
                   + 2 * K_ym_b * u[-2, 1:-1]) / (2.0 * dx2)
                - denom_b * u[-1, 1:-1])))
            worst = max(worst, np.max(np.abs(
                num_const_l
                + (K_yp_l * u[2:, 0] + K_ym_l * u[:-2, 0]
                   + 2 * K_xp_l * u[1:-1, 1]) / (2.0 * dx2)
                - denom_l * u[1:-1, 0])))
            worst = max(worst, np.max(np.abs(
                num_const_r
                + (K_yp_r * u[2:, -1] + K_ym_r * u[:-2, -1]
                   + 2 * K_xm_r * u[1:-1, -2]) / (2.0 * dx2)
                - denom_r * u[1:-1, -1])))
            for (y, x, ny_nb, nx_nb, Ky, Kx, denom_c, num_c_const) in corner_data:
                worst = max(worst, abs(
                    num_c_const
                    + (2 * Ky * u[ny_nb, x] + 2 * Kx * u[y, nx_nb]) / (2.0 * dx2)
                    - denom_c * u[y, x]))
            return float(worst) / rhs_scale

        iteration = 0
        converged = False
        guard_attempts = 0
        residual = float("inf")

        while iteration < max_iterations:
            # FIX C-5: cooperative cancellation, polled every iteration so the
            # Stop button and window close take effect immediately.
            if should_continue is not None and not should_continue():
                logger.info(
                    "Steady-state solve interrupted after %d iteration(s) "
                    "(relative residual %.3e).", iteration, residual,
                )
                self.u = u
                return self.u, iteration, False

            # One red sweep then one black sweep, boundaries included (FIX H-1).
            for blocks in colour_blocks:
                for block in blocks:
                    block(omega)

            iteration += 1

            if iteration % check_interval and iteration != max_iterations:
                continue

            diverged = not np.isfinite(u).all()
            if not diverged:
                residual = relative_residual()
                diverged = not np.isfinite(residual) or residual > 1e12

            if diverged:
                # FIX C-4: bounded number of restarts, then give up honestly.
                guard_attempts += 1
                if guard_attempts > max_guard_attempts:
                    logger.error(
                        "Steady-state solve diverged again after %d guard "
                        "restart(s); giving up. Check the k and h fields for "
                        "extreme values.", max_guard_attempts,
                    )
                    self.u = u
                    return self.u, iteration, False
                logger.warning(
                    "Steady-state solve diverged at omega=%.2f (attempt %d of "
                    "%d); restarting from ambient with omega=1.0.",
                    omega, guard_attempts, max_guard_attempts,
                )
                u[:] = self.T_amb
                omega = 1.0
                iteration = 0
                residual = float("inf")
                continue

            if residual <= tolerance:
                converged = True
                break

        if converged:
            logger.info(
                "Steady-state converged in %d iteration(s); relative residual "
                "%.3e.", iteration, residual,
            )
        else:
            logger.warning(
                "Steady-state did NOT converge: stopped at the %d-iteration "
                "cap with relative residual %.3e (target %.1e). The reported "
                "temperatures are not a solution - treat them as a lower "
                "bound.", iteration, residual, tolerance,
            )

        self.u = u
        return self.u, iteration, converged
