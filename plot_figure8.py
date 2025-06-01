from instability_curves import InstabilityCurves
import numpy as np
import matplotlib.pyplot as plt

state = 'B-'
eta_val = 2.0139
k_val = 2
q_val = 0.1
Nt = 4000

model = InstabilityCurves(state=state, eta_val=eta_val, k_val=k_val, q_val=q_val, Nt=Nt)

eps2_bound = model.generate_epsilon2_values()
n_rows = len(eps2_bound)
n_cols = 3  

fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3))

if n_rows == 1:
    axes = axes.reshape(1, 3)

for i, epsilon_ in enumerate(eps2_bound):
    a, b = epsilon_[0], epsilon_[1]
    epsilon2 = 0.5 * (a + b)
    print(f"Running simulation for {a:.4f} < ε² = {epsilon2:.4f} < {b:.4f}")

    model.run_simulation(epsilon2)

    A = model.A_result
    B = model.B_result
    F_shifted = model.F_shifted

    axes[i, 0].imshow(A, cmap='viridis', origin='lower')
    axes[i, 0].set_title("A", fontsize=9)
    axes[i, 0].axis('off')

    axes[i, 1].imshow(B, cmap='viridis', origin='lower')
    axes[i, 1].set_title("B", fontsize=9)
    axes[i, 1].axis('off')

    axes[i, 2].imshow(np.log1p(F_shifted), cmap='gray')
    axes[i, 2].set_title("FFT log(|B - B⁻|)", fontsize=9)
    axes[i, 2].axis('off')

    axes[i, 0].text(-0.1, 0.5, f"{a:.2f} < ε² = {epsilon2:.2f} < {b:.2f}",
                    transform=axes[i, 0].transAxes, va='center', ha='right', fontsize=8)

fig.suptitle(
    f"State = {state}, η = {eta_val}, k = {k_val}, q = {q_val}, Nt = {Nt} \n"
    f"lower bound given by paper = {eps2_bound[0][1]} \n"
    f"lower bound given by us = {eps2_bound[1][1]} \n"
    f"upper bound given by paper = {eps2_bound[2][1]}",
    fontsize=12, y=0.99
)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
