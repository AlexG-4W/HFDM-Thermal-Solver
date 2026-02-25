import numpy as np
import config

class PCBSolver:
    def __init__(self, nx, ny, Q_matrix, material_name=None, h_matrix=None, layers=None, copper_oz=None, K_matrix=None):
        self.nx = nx
        self.ny = ny
        self.Q = Q_matrix
        
        # Initialize temperature field with ambient temperature
        self.u = np.full((ny, nx), config.T_amb, dtype=np.float64)
        
        # Determine base properties
        if material_name:
            props = config.MATERIALS.get(material_name)
            if not props:
                raise ValueError(f"Material {material_name} not found.")
            self.rho = props["rho"]
            self.cp = props["cp"]
            base_k = props["k"]
        else:
            self.rho = config.rho
            self.cp = config.cp
            base_k = config.MATERIALS.get(config.SELECTED_MATERIAL, config.MATERIALS["FR-4"])["k"]

        # If K_matrix is provided, use it. Otherwise, use k_eff logic.
        if K_matrix is not None:
            if K_matrix.shape != (ny, nx):
                raise ValueError(f"K_matrix shape {K_matrix.shape} must match grid {(ny, nx)}")
            self.K = K_matrix
            self.k_eff = np.mean(K_matrix)
        else:
            l = layers if layers is not None else config.BOARD_LAYERS
            oz = copper_oz if copper_oz is not None else config.COPPER_OZ
            k = config.calculate_k_eff(layers=l, copper_oz=oz, substrate_k=base_k)
            self.k_eff = k
            self.K = np.full((ny, nx), k)

        # Stability based on max conductivity
        k_max = np.max(self.K)
        alpha_max = k_max / (self.rho * self.cp)
        self.dt = 0.9 * (config.dx**2) / (4 * alpha_max)
            
        # Convective heat transfer coefficient (h) - Global or Matrix
        if h_matrix is not None:
            if h_matrix.shape != (ny, nx):
                raise ValueError(f"h_matrix shape {h_matrix.shape} must match grid {(ny, nx)}")
            self.h = h_matrix
        else:
            self.h = np.full((ny, nx), config.h)
            
        # Simulation constants
        self.coeff = self.dt / (self.rho * self.cp * config.dx**2)
        self.Source_Term = (self.Q * self.dt) / (self.rho * self.cp)
        self.Surface_Loss_Coeff = (2 * self.h * self.dt) / (self.rho * self.cp * config.d)
        
        # Precompute interface conductivities for interior nodes (1:-1, 1:-1)
        # These are (K_i + K_neighbor) / 2.0
        self.K_xp = (self.K[1:-1, 1:-1] + self.K[1:-1, 2:]) / 2.0
        self.K_xm = (self.K[1:-1, 1:-1] + self.K[1:-1, :-2]) / 2.0
        self.K_yp = (self.K[1:-1, 1:-1] + self.K[2:, 1:-1]) / 2.0
        self.K_ym = (self.K[1:-1, 1:-1] + self.K[:-2, 1:-1]) / 2.0

    def update_q_matrix(self, new_Q):
        """Updates the heat generation matrix and source term."""
        self.Q = new_Q
        self.Source_Term = (self.Q * self.dt) / (self.rho * self.cp)

    def step(self):
        """Performs one time step using a fully vectorized explicit Euler scheme for heterogeneous K."""
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

        # 3. Boundary Edges (Convection at board edges)
        # Top Edge (y=0)
        K_top = self.K[0, 1:-1]
        K_xp_t = (K_top + self.K[0, 2:]) / 2.0
        K_xm_t = (K_top + self.K[0, :-2]) / 2.0
        K_yp_t = (K_top + self.K[1, 1:-1]) / 2.0
        flux_x_t = K_xp_t * (u[0, 2:] - u[0, 1:-1]) + K_xm_t * (u[0, :-2] - u[0, 1:-1])
        flux_y_t = 2 * K_yp_t * (u[1, 1:-1] - u[0, 1:-1]) + 2 * self.h[0, 1:-1] * dx * (T_amb - u[0, 1:-1])
        u_new[0, 1:-1] += self.coeff * (flux_x_t + flux_y_t)

        # Bottom Edge (y=ny-1)
        K_bot = self.K[-1, 1:-1]
        K_xp_b = (K_bot + self.K[-1, 2:]) / 2.0
        K_xm_b = (K_bot + self.K[-1, :-2]) / 2.0
        K_ym_b = (K_bot + self.K[-2, 1:-1]) / 2.0
        flux_x_b = K_xp_b * (u[-1, 2:] - u[-1, 1:-1]) + K_xm_b * (u[-1, :-2] - u[-1, 1:-1])
        flux_y_b = 2 * K_ym_b * (u[-2, 1:-1] - u[-1, 1:-1]) + 2 * self.h[-1, 1:-1] * dx * (T_amb - u[-1, 1:-1])
        u_new[-1, 1:-1] += self.coeff * (flux_x_b + flux_y_b)

        # Left Edge (x=0)
        K_left = self.K[1:-1, 0]
        K_yp_l = (K_left + self.K[2:, 0]) / 2.0
        K_ym_l = (K_left + self.K[:-2, 0]) / 2.0
        K_xp_l = (K_left + self.K[1:-1, 1]) / 2.0
        flux_y_l = K_yp_l * (u[2:, 0] - u[1:-1, 0]) + K_ym_l * (u[:-2, 0] - u[1:-1, 0])
        flux_x_l = 2 * K_xp_l * (u[1:-1, 1] - u[1:-1, 0]) + 2 * self.h[1:-1, 0] * dx * (T_amb - u[1:-1, 0])
        u_new[1:-1, 0] += self.coeff * (flux_y_l + flux_x_l)

        # Right Edge (x=nx-1)
        K_right = self.K[1:-1, -1]
        K_yp_r = (K_right + self.K[2:, -1]) / 2.0
        K_ym_r = (K_right + self.K[:-2, -1]) / 2.0
        K_xm_r = (K_right + self.K[1:-1, -2]) / 2.0
        flux_y_r = K_yp_r * (u[2:, -1] - u[1:-1, -1]) + K_ym_r * (u[:-2, -1] - u[1:-1, -1])
        flux_x_r = 2 * K_xm_r * (u[1:-1, -2] - u[1:-1, -1]) + 2 * self.h[1:-1, -1] * dx * (T_amb - u[1:-1, -1])
        u_new[1:-1, -1] += self.coeff * (flux_y_r + flux_x_r)

        # 4. Corners
        corners = [(0, 0, 1, 1), (0, -1, 1, -2), (-1, 0, -2, 1), (-1, -1, -2, -2)]
        for (y, x, ny_nb, nx_nb) in corners:
            K_c = self.K[y, x]
            K_nb_y = (K_c + self.K[ny_nb, x]) / 2.0
            K_nb_x = (K_c + self.K[y, nx_nb]) / 2.0
            flux_c = 2 * K_nb_y * (u[ny_nb, x] - u[y, x]) + 2 * K_nb_x * (u[y, nx_nb] - u[y, x]) + 4 * self.h[y, x] * dx * (T_amb - u[y, x])
            u_new[y, x] += self.coeff * flux_c

        self.u = u_new
        return self.u

    def get_probe_temperatures(self, indices):
        """Returns the temperatures at the specified (y, x) indices."""
        return [self.u[y, x] for y, x in indices]

    def solve_steady_state(self, tolerance=1e-4, max_iterations=50000):
        """Finds the final temperature distribution using Jacobi iteration with heterogeneous K."""
        dx = config.dx
        dx2 = dx**2
        d = config.d
        h = self.h
        T_amb = config.T_amb
        K = self.K
        
        # Precompute interface conductivities (Summed form for Jacobi)
        # Interior: (K_xp + K_xm + K_yp + K_ym)
        K_xp = (K[1:-1, 1:-1] + K[1:-1, 2:])
        K_xm = (K[1:-1, 1:-1] + K[1:-1, :-2])
        K_yp = (K[1:-1, 1:-1] + K[2:, 1:-1])
        K_ym = (K[1:-1, 1:-1] + K[:-2, 1:-1])
        
        # Denominator for interior: Sum(K_interface)/2dx2 + 2h/d
        denom_int = (K_xp + K_xm + K_yp + K_ym) / (2.0 * dx2) + (2 * h[1:-1, 1:-1] / d)
        num_const_int = self.Q[1:-1, 1:-1] + (2 * h[1:-1, 1:-1] * T_amb / d)
        
        iteration = 0
        error = 1.0
        
        while error > tolerance and iteration < max_iterations:
            u = self.u
            u_new = np.copy(u)
            
            # 1. Interior Nodes Jacobi Update
            S_int = (K_xp * u[1:-1, 2:] + K_xm * u[1:-1, :-2] + K_yp * u[2:, 1:-1] + K_ym * u[:-2, 1:-1]) / (2.0 * dx2)
            u_new[1:-1, 1:-1] = (S_int + num_const_int) / denom_int
            
            # 2. Boundary Edges (Ghost Node method adapted for steady state with variable K)
            # Top Edge (y=0)
            K_t = K[0, 1:-1]
            K_xp_t = (K_t + K[0, 2:])
            K_xm_t = (K_t + K[0, :-2])
            K_yp_t = (K_t + K[1, 1:-1])
            denom_t = (K_xp_t + K_xm_t + 2 * K_yp_t) / (2.0 * dx2) + (2 * h[0, 1:-1] / d) + (2 * h[0, 1:-1] / dx)
            num_t = (K_xp_t * u[0, 2:] + K_xm_t * u[0, :-2] + 2 * K_yp_t * u[1, 1:-1]) / (2.0 * dx2) + \
                    self.Q[0, 1:-1] + (2 * h[0, 1:-1] * T_amb / d) + (2 * h[0, 1:-1] * T_amb / dx)
            u_new[0, 1:-1] = num_t / denom_t

            # Bottom Edge (y=ny-1)
            K_b = K[-1, 1:-1]
            K_xp_b = (K_b + K[-1, 2:])
            K_xm_b = (K_b + K[-1, :-2])
            K_ym_b = (K_b + K[-2, 1:-1])
            denom_b = (K_xp_b + K_xm_b + 2 * K_ym_b) / (2.0 * dx2) + (2 * h[-1, 1:-1] / d) + (2 * h[-1, 1:-1] / dx)
            num_b = (K_xp_b * u[-1, 2:] + K_xm_b * u[-1, :-2] + 2 * K_ym_b * u[-2, 1:-1]) / (2.0 * dx2) + \
                    self.Q[-1, 1:-1] + (2 * h[-1, 1:-1] * T_amb / d) + (2 * h[-1, 1:-1] * T_amb / dx)
            u_new[-1, 1:-1] = num_b / denom_b

            # Left Edge (x=0)
            K_l = K[1:-1, 0]
            K_yp_l = (K_l + K[2:, 0])
            K_ym_l = (K_l + K[:-2, 0])
            K_xp_l = (K_l + K[1:-1, 1])
            denom_l = (K_yp_l + K_ym_l + 2 * K_xp_l) / (2.0 * dx2) + (2 * h[1:-1, 0] / d) + (2 * h[1:-1, 0] / dx)
            num_l = (K_yp_l * u[2:, 0] + K_ym_l * u[:-2, 0] + 2 * K_xp_l * u[1:-1, 1]) / (2.0 * dx2) + \
                    self.Q[1:-1, 0] + (2 * h[1:-1, 0] * T_amb / d) + (2 * h[1:-1, 0] * T_amb / dx)
            u_new[1:-1, 0] = num_l / denom_l

            # Right Edge (x=nx-1)
            K_r = K[1:-1, -1]
            K_yp_r = (K_r + K[2:, -1])
            K_ym_r = (K_r + K[:-2, -1])
            K_xm_r = (K_r + K[1:-1, -2])
            denom_r = (K_yp_r + K_ym_r + 2 * K_xm_r) / (2.0 * dx2) + (2 * h[1:-1, -1] / d) + (2 * h[1:-1, -1] / dx)
            num_r = (K_yp_r * u[2:, -1] + K_ym_r * u[:-2, -1] + 2 * K_xm_r * u[1:-1, -2]) / (2.0 * dx2) + \
                    self.Q[1:-1, -1] + (2 * h[1:-1, -1] * T_amb / d) + (2 * h[1:-1, -1] * T_amb / dx)
            u_new[1:-1, -1] = num_r / denom_r

            # 3. Corners (Double Ghost Nodes)
            corners = [(0, 0, 1, 1), (0, -1, 1, -2), (-1, 0, -2, 1), (-1, -1, -2, -2)]
            for (y, x, ny_nb, nx_nb) in corners:
                K_c = K[y, x]
                K_nb_y = (K_c + K[ny_nb, x])
                K_nb_x = (K_c + K[y, nx_nb])
                denom_c = (2 * K_nb_y + 2 * K_nb_x) / (2.0 * dx2) + (2 * h[y, x] / d) + (4 * h[y, x] / dx)
                num_c = (2 * K_nb_y * u[ny_nb, x] + 2 * K_nb_x * u[y, nx_nb]) / (2.0 * dx2) + \
                        self.Q[y, x] + (2 * h[y, x] * T_amb / d) + (4 * h[y, x] * T_amb / dx)
                u_new[y, x] = num_c / denom_c

            error = np.max(np.abs(u_new - u))
            self.u = u_new
            iteration += 1
            
        return self.u, iteration
