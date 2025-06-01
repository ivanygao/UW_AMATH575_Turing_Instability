from instability_curves import InstabilityCurves
import numpy as np
import matplotlib.pyplot as plt

model = InstabilityCurves(state='B-', k_val=2, q_val=0)

eta_vals = np.linspace(1.8, 2.3, 300)

tau_determinant = [model.solve(eta, mode='tau1') for eta in eta_vals]
tau_trace = [model.solve(eta, mode='tau2') for eta in eta_vals]

plt.figure(figsize=(10, 6))
plt.plot(eta_vals, tau_determinant, label='Turing curve', color='red', linewidth=2)
plt.plot(eta_vals, tau_trace, label='Hopf curve', color='black', linestyle='--', linewidth=2)

plt.fill_between(eta_vals, 1, 4, where=(eta_vals < 2), color='pink', alpha=0.3)
plt.text(1.85, 3.5, 'SRP', color='red', fontsize=12, rotation=45)
plt.text(2.05, 2.5, 'Turing space', color='black', fontsize=12, rotation=45)
plt.text(1.85, 1.5, 'Trivial state (1,0)', color='black', fontsize=10, rotation=90)
plt.text(2.2, 1.5, f'State ({model.state})', color='black', fontsize=10)

plt.xlabel('η')
plt.ylabel('τ')
plt.title('Turing and Hopf Instability Curves')
plt.grid(alpha=0.3)
plt.legend()
plt.ylim(1, 4)
plt.xlim(min(eta_vals), max(eta_vals))
plt.tight_layout()
plt.show()
