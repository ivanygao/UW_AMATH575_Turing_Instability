import numpy as np
import matplotlib.pyplot as plt
from instability_curves import InstabilityCurves

model = InstabilityCurves(state='B+') 

eta = 2.01
τ = 4.99
epsilon_values = [3.0, 4.1, 5.0]
k_values = [e / τ for e in epsilon_values]
q_vals = np.linspace(0, 4, 500)
q2_vals = q_vals**2


plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
model.state = 'B+'
model._define_symbols()
model._compute_symbolic_expressions()
model._prepare_numeric_functions()

for epsilon_sq in epsilon_values:
    k = epsilon_sq/τ
    for q in np.linspace(0, 4, 500):
        y_vals = [model._f_det_fixed_eta(eta, q2, τ, k) for q2 in q2_vals]
    color = 'black' if epsilon_sq == 3.0 else 'red' if epsilon_sq == 4.1 else 'green'
    plt.plot(q2_vals, y_vals, label=f'ε² = {epsilon_sq}', color=color)

plt.axhline(0, color='black', linestyle='--')
plt.text(0.2, 5.5, 'η = 2.01', fontsize=10)
plt.text(0.2, 5.0, 'τ = 4.99', fontsize=10)
plt.title('(a) B⁺ state')
plt.xlabel('$q^2$')
plt.ylabel('$|M(q^2)|$')
plt.legend()
plt.grid()
plt.ylim(-4, 6)
plt.xlim(0, 4)

plt.subplot(1, 2, 2)
model.state = 'B-'
model._define_symbols()
model._compute_symbolic_expressions()
model._prepare_numeric_functions()

for epsilon_sq in epsilon_values:
    k = epsilon_sq/τ
    for q in np.linspace(0, 4, 500):
        y_vals = [model._f_det_fixed_eta(eta, q2, τ, k) for q2 in q2_vals]
    color = 'black' if epsilon_sq == 3.0 else 'red' if epsilon_sq == 4.1 else 'green'
    plt.plot(q2_vals, y_vals, label=f'ε² = {epsilon_sq}', color=color)


plt.axhline(0, color='black', linestyle='--')
plt.fill_between(q2_vals, -4, 0, where=(q2_vals >= 1.2) & (q2_vals <= 2.0), color='grey', alpha=0.3)
plt.vlines(1.2, -4, 0, colors='black', linestyles='dashed')
plt.vlines(2.0, -4, 0, colors='black', linestyles='dashed')
plt.text(1.5, -3.5, 'Patterns', fontsize=10, color='green')
plt.text(2.5, 4.0, 'No patterns', fontsize=10, color='red', fontstyle='italic')
plt.text(1.7, 0.2, 'Limit of stability', fontsize=10, color='red', fontstyle='italic')
plt.text(0.2, 5.5, 'η = 2.01', fontsize=10)
plt.text(0.2, 5.0, 'τ = 4.99', fontsize=10)
plt.title('(b) B⁻ state')
plt.xlabel('$q^2$')
plt.ylabel('$|M(q^2)|$')
plt.legend()
plt.grid()
plt.ylim(-4, 6)
plt.xlim(0, 4)

plt.tight_layout()
plt.show()


