import numpy as np
import config


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
    def __init__(self, nx, ny, Q_matrix, material_name=None, h_matrix=None,
                 layers=None, copper_oz=None, K_matrix=None):
        self.nx = nx
        self.ny = ny
        self.Q = Q_matrix

        # Initialize temperature field with ambient temperature
        self.u = np.full((ny, nx), config.T_amb, dtype=np.float64)

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
            self.h = h_matrix
        else:
            self.h = np.full((ny, nx), config.h)

        # FIX S-6: Full CFL condition including convection surface sink.
        # dt < 1 / (4*alpha_max/dx^2 + 2*h_max/(rho*cp*d))
        k_max = np.max(self.K)
        alpha_max = k_max / (self.rho * self.cp)
        h_max = np.max(self.h)
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
        T_amb = config.T_amb

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
        return self.u

    def get_probe_temperatures(self, indices):
        """Returns the temperatures at the specified (y, x) indices."""
        return [self.u[y, x] for y, x in indices]

    def solve_steady_state(self, tolerance=1e-4, max_iterations=50000, omega=1.5):
        """
        Finds the final temperature distribution using Red-Black SOR
        with heterogeneous K and harmonic-mean interface conductivities.

        FIX S-1: Harmonic mean for all interface K.
        FIX S-3: Red-Black SOR replaces Jacobi (~5-10x faster convergence).
        FIX S-5: Pre-allocated buffers, no np.copy per iteration.

        Parameters
        ----------
        omega : float
            SOR relaxation factor. Values 1.0 < omega < 2.0.
            omega=1.0 is equivalent to Gauss-Seidel.
            Default 1.5 is safe for heterogeneous h fields.
            If divergence is detected, automatically falls back to omega=1.0.
        """
        dx = config.dx
        dx2 = dx**2
        d = config.d
        h = self.h
        T_amb = config.T_amb
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

        # FIX S-3: Red-Black SOR iteration masks.
        # "Red" cells: (i+j) % 2 == 0, "Black" cells: (i+j) % 2 == 1
        # For the interior sub-grid [1:-1, 1:-1], indices run from 0..ny-3, 0..nx-3
        ny_int = self.ny - 2
        nx_int = self.nx - 2
        ii, jj = np.mgrid[0:ny_int, 0:nx_int]
        # Offset by (1,1) to get global parity: global_i = ii+1, global_j = jj+1
        red_mask = ((ii + jj) % 2 == 0)
        black_mask = ~red_mask

        iteration = 0
        error = 1.0

        # FIX S-5: Pre-allocate error buffer
        u_prev = np.empty_like(u)

        while error > tolerance and iteration < max_iterations:
            u_prev[:] = u

            # Divergence guard: if NaN appears, restart with omega=1.0
            if iteration > 0 and (np.any(np.isnan(u)) or error > 1e10):
                u[:] = config.T_amb
                u_prev[:] = u
                omega = 1.0
                iteration = 0
                error = 1.0
                continue

            # --- Red sweep (interior) ---
            S_red = (K_xp * u[1:-1, 2:] + K_xm * u[1:-1, :-2]
                     + K_yp * u[2:, 1:-1] + K_ym * u[:-2, 1:-1]) / (2.0 * dx2)
            u_gs = (S_red + num_const_int) / denom_int
            u[1:-1, 1:-1] = np.where(red_mask,
                                      u[1:-1, 1:-1] + omega * (u_gs - u[1:-1, 1:-1]),
                                      u[1:-1, 1:-1])

            # --- Black sweep (interior) ---
            S_black = (K_xp * u[1:-1, 2:] + K_xm * u[1:-1, :-2]
                       + K_yp * u[2:, 1:-1] + K_ym * u[:-2, 1:-1]) / (2.0 * dx2)
            u_gs = (S_black + num_const_int) / denom_int
            u[1:-1, 1:-1] = np.where(black_mask,
                                      u[1:-1, 1:-1] + omega * (u_gs - u[1:-1, 1:-1]),
                                      u[1:-1, 1:-1])

            # --- Boundary Edges (SOR applied uniformly) ---
            # Top Edge
            num_t = (K_xp_t * u[0, 2:] + K_xm_t * u[0, :-2] + 2 * K_yp_t * u[1, 1:-1]) / (2.0 * dx2) + num_const_t
            u_gs_t = num_t / denom_t
            u[0, 1:-1] += omega * (u_gs_t - u[0, 1:-1])

            # Bottom Edge
            num_b = (K_xp_b * u[-1, 2:] + K_xm_b * u[-1, :-2] + 2 * K_ym_b * u[-2, 1:-1]) / (2.0 * dx2) + num_const_b
            u_gs_b = num_b / denom_b
            u[-1, 1:-1] += omega * (u_gs_b - u[-1, 1:-1])

            # Left Edge
            num_l = (K_yp_l * u[2:, 0] + K_ym_l * u[:-2, 0] + 2 * K_xp_l * u[1:-1, 1]) / (2.0 * dx2) + num_const_l
            u_gs_l = num_l / denom_l
            u[1:-1, 0] += omega * (u_gs_l - u[1:-1, 0])

            # Right Edge
            num_r = (K_yp_r * u[2:, -1] + K_ym_r * u[:-2, -1] + 2 * K_xm_r * u[1:-1, -2]) / (2.0 * dx2) + num_const_r
            u_gs_r = num_r / denom_r
            u[1:-1, -1] += omega * (u_gs_r - u[1:-1, -1])

            # Corners
            for (y, x, ny_nb, nx_nb, Ky, Kx, denom_c, num_c_const) in corner_data:
                num_c = (2 * Ky * u[ny_nb, x] + 2 * Kx * u[y, nx_nb]) / (2.0 * dx2) + num_c_const
                u_gs_c = num_c / denom_c
                u[y, x] += omega * (u_gs_c - u[y, x])

            error = np.max(np.abs(u - u_prev))
            iteration += 1

        self.u = u
        return self.u, iteration
