import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftshift

class InstabilityCurves:
    def __init__(self, state='B-', eta_val=2.0139, k_val= 2, q_val=0, tau_val=2.7339,dt=0.01, Nt=40000, size=200, phi = 0.01, qx = 1, qy = 1):
        self.k_val = k_val
        self.state = state
        self.q_val = q_val
        self.eta_vals = eta_val
        self.tau_val = tau_val
        self.dt = dt
        self.Nt = Nt
        self.size = size
        self.phi = phi
        self.qx = qx
        self.qy = qy

        self.eta = sp.Symbol('eta', real=True, positive=True)
        self.tau = sp.Symbol('tau', real=True, positive=True)
        self.q_sq = sp.Symbol('q_sq', real=True, positive=True)
        self.k = sp.Symbol('k', real=True, positive=True)
        self.epsilon_sq = self.k * self.tau
        self.sqrt_term = sp.sqrt(self.eta**2 - 4)

        self._define_symbols()
        self._compute_symbolic_expressions()
        self._prepare_numeric_functions()

    def _define_symbols(self):
        if self.state == 'B+':
            self.A_e = (self.eta + self.sqrt_term) / (2 * self.eta)
            self.B_e = (self.eta - self.sqrt_term) / 2
        elif self.state == 'B-':
            self.A_e = (self.eta - self.sqrt_term) / (2 * self.eta)
            self.B_e = (self.eta + self.sqrt_term) / 2
        else:
            raise ValueError("State must be 'B+' or 'B-'")

    def _compute_symbolic_matrix(self):
        A_e, B_e = self.A_e, self.B_e
        eta, tau, q_sq = self.eta, self.tau, self.q_sq
        epsilon_sq = self.epsilon_sq

        M11 = tau * (2 * eta * A_e * B_e - 1 - q_sq / epsilon_sq)
        M12 = tau * eta * B_e**2
        M21 = -2 * A_e * B_e
        M22 = -(q_sq + B_e**2 + 1)

        M = sp.Matrix([[M11, M12], [M21, M22]])

        self.M_sym = M
        self.trace_expr = M.trace()
        self.det_expr = M.det()

    def _compute_symbolic_roots(self):
        tau, q_sq = self.tau, self.q_sq
        self.q_sq_det_expr = sp.solve(self.det_expr, q_sq)[0]
        self.tau_trace_expr = sp.solve(self.trace_expr, tau)[0]

    def _compute_symbolic_discriminant(self):
        tau, q_sq = self.tau, self.q_sq
        poly = sp.Poly(self.det_expr, q_sq)
        a, b, c = poly.all_coeffs()
        self.delta = b**2 - 4 * a * c
        self.tau_delta_expr = sp.solve(self.delta , tau)[1]

    def _compute_symbolic_expressions(self):
        self._compute_symbolic_matrix()
        self._compute_symbolic_roots()
        self._compute_symbolic_discriminant()

    def _prepare_numeric_functions(self):
        self.f_tau_trace = sp.lambdify((self.eta, self.q_sq, self.k), self.tau_trace_expr, modules='numpy')
        self.f_tau_delta = sp.lambdify((self.eta, self.q_sq, self.k), self.tau_delta_expr, modules='numpy')
        self.f_A = sp.lambdify(self.eta, self.A_e, modules='numpy')
        self.f_B = sp.lambdify(self.eta, self.B_e, modules='numpy')
        self._f_det_fixed_eta = sp.lambdify((self.eta, self.q_sq, self.tau, self.k), self.det_expr, modules='numpy')

    def solve(self, eta_val, mode='tau1'):
        if eta_val**2 < 4:
            return np.nan
        q_sq_val = self.q_val**2
        if mode == 'tau1':
            return self.f_tau_delta(eta_val, q_sq_val, self.k_val)
        elif mode == 'tau2':
            return self.f_tau_trace(eta_val, q_sq_val, self.k_val)
        else:
            raise ValueError("mode must be 'tau1' or 'tau2'")

    @staticmethod
    def laplacian(Z):
            return (np.roll(Z, 1, axis=0) + np.roll(Z, -1, axis=0) +
                    np.roll(Z, 1, axis=1) + np.roll(Z, -1, axis=1) - 4 * Z)

    def run_simulation(self, epsilon2):
        eta = self.eta_vals
        tau = self.tau_val
        dt = self.dt
        Nt = self.Nt
        size = self.size
        phi = self.phi
        qx = self.qx
        qy = self.qy

        A_minus = self.f_A(eta)
        B_minus = self.f_B(eta)

        x = np.linspace(0, size, size, endpoint=False)
        y = np.linspace(0, size, size, endpoint=False)
        X, Y = np.meshgrid(x, y, indexing='ij')
        perturb = phi * np.cos(qx * X) * np.cos(qy * Y)

        A = A_minus + perturb
        B = B_minus + perturb


        for n in range(Nt):
            LA = self.laplacian(A)
            LB = self.laplacian(B)
            A += (LA + (1 - A) - B**2 * A) * dt
            B += tau * (eta * B**2 * A - B + (1 / epsilon2) * LB) * dt

        self.A_result = A
        self.B_result = B

        self.F = fft2(B - B_minus)
        self.F_shifted = fftshift(np.abs(self.F))

    def batch_simulate(self, epsilon2_list):
        for eps in epsilon2_list:
            print(f"Running ε² = {eps}")
            self.run_simulation(eps)

    def generate_epsilon2_values(self):
        eta = self.eta_vals
        k = self.k_val

        B_e = self.f_B(eta)
        tau_crit = self.f_tau_delta(eta, self.q_val**2, k)
        eps2_tau_crit = k * tau_crit
        eps2_be_eta = B_e * eta
        eps2_empirical = 11.65 * B_e * eta - 8

        eps2_key = [
            0,
            float(eps2_be_eta),
            float(eps2_tau_crit),
            float(eps2_empirical),
            100
        ]

        eps2_key = sorted(set(eps2_key))

        eps2_bound = []
        for i in range(len(eps2_key) - 1):
            eps2_bound.append([eps2_key[i], eps2_key[i + 1]]) 

        return eps2_bound
        