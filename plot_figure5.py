from instability_curves import InstabilityCurves
import numpy as np
import matplotlib.pyplot as plt

eta_values = np.linspace(0, 15, 300)
A_stable_vals, B_stable_vals, A_unstable_vals, B_unstable_vals = [], [], [], []

sim_eta = np.linspace(2, 14, 50)
sim_model = InstabilityCurves(state='B-', q_val=0) 
sim_A = sim_model.f_A(sim_eta)
sim_B = sim_model.f_B(sim_eta)

model_stable = InstabilityCurves(state='B-', q_val=0)
model_unstable = InstabilityCurves(state='B+',  q_val=0)

for eta in eta_values:
    if eta >= 2:
        A_stable_vals.append(model_stable.f_A(eta))
        B_stable_vals.append(model_stable.f_B(eta))
        A_unstable_vals.append(model_unstable.f_A(eta))
        B_unstable_vals.append(model_unstable.f_B(eta))
    else:
        A_stable_vals.append(1.0)
        B_stable_vals.append(0.0)
        A_unstable_vals.append(np.nan)
        B_unstable_vals.append(np.nan)

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(eta_values, A_stable_vals, label='Stable A', color='black')
plt.plot(eta_values, A_unstable_vals, '--', label='Unstable A', color='black')
plt.scatter(sim_eta, sim_A, marker='*', color='black', label='Simulations')
plt.title('Figure 5 (a) - A vs η')
plt.xlabel('η')
plt.ylabel('A')
plt.legend()
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(eta_values, B_stable_vals, label='Stable B', color='black')
plt.plot(eta_values, B_unstable_vals, '--', label='Unstable B', color='black')
plt.scatter(sim_eta, sim_B, marker='*', color='black', label='Simulations')
plt.title('Figure 5 (b) - B vs η')
plt.xlabel('η')
plt.ylabel('B')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
