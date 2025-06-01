import numpy as np
import matplotlib.pyplot as plt
from instability_curves import InstabilityCurves  # 假设你将类保存为 instability_curves.py

# 初始化类
model = InstabilityCurves(
    state='B-', eta_val=2.0139, k_val=2, q_val=0,
    tau_val=2.7339, dt=0.01, Nt=40000, size=200,
    phi=0.01, qx=1, qy=1
)

# === LBM 模拟函数 ===
def run_lbm_simulation(model, epsilon2):
    size = model.size
    Nt = model.Nt
    dt = model.dt
    eta = model.eta_vals
    tau = model.tau_val
    phi = model.phi
    qx = model.qx
    qy = model.qy

    A0 = model.f_A(eta)
    B0 = model.f_B(eta)

    # D2Q9 格式定义
    c = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
                  [1, 1], [-1, 1], [-1, -1], [1, -1]])
    w = np.array([4/9] + [1/9]*4 + [1/36]*4)

    # 初始化
    A = np.full((size, size), A0)
    B = np.full((size, size), B0)
    X, Y = np.meshgrid(np.arange(size), np.arange(size), indexing='ij')
    perturb = phi * np.cos(qx * X * 2 * np.pi / size) * np.cos(qy * Y * 2 * np.pi / size)
    A += perturb
    B += perturb

    fA = np.array([w[i] * A for i in range(9)])
    fB = np.array([w[i] * B for i in range(9)])

    tauA = 3 * 1 + 0.5   # D_A = 1 in lattice units
    tauB = 3 * (1 / epsilon2) + 0.5

    for t in range(Nt):
        rhoA = np.sum(fA, axis=0)
        rhoB = np.sum(fB, axis=0)

        # 反应项
        RA = -rhoA * rhoB**2 + (1 - rhoA)
        RB = eta * rhoB**2 * rhoA - rhoB

        # 平衡态
        feqA = np.array([w[i] * rhoA for i in range(9)])
        feqB = np.array([w[i] * rhoB for i in range(9)])

        fA += -(fA - feqA) / tauA + dt * np.array([w[i] * RA for i in range(9)])
        fB += -(fB - feqB) / tauB + dt * np.array([w[i] * RB for i in range(9)])

        for i in range(9):
            fA[i] = np.roll(fA[i], c[i][0], axis=1)
            fA[i] = np.roll(fA[i], c[i][1], axis=0)
            fB[i] = np.roll(fB[i], c[i][0], axis=1)
            fB[i] = np.roll(fB[i], c[i][1], axis=0)

    model.B_result_lbm = np.sum(fB, axis=0)
    model.F_lbm = np.fft.fft2(model.B_result_lbm - B0)
    model.F_shifted_lbm = np.fft.fftshift(np.abs(model.F_lbm))

# === 主程序调用 ===
epsilon2_ranges = model.generate_epsilon2_values()

for low, high in epsilon2_ranges:
    eps2 = 0.5 * (low + high)
    print(f"Running LBM for ε² = {eps2:.4f}")
    run_lbm_simulation(model, eps2)

    # 绘制图像
    plt.figure(figsize=(6, 5))
    plt.imshow(model.B_result_lbm, cmap='inferno', origin='lower')
    plt.colorbar(label='B Concentration')
    plt.title(f'LBM: Pattern at $\epsilon^2$ = {eps2:.4f}')
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(6, 5))
    plt.imshow(np.log1p(model.F_shifted_lbm), cmap='magma', origin='lower')
    plt.colorbar(label='log(|F|)')
    plt.title(f'LBM: Fourier Spectrum at $\epsilon^2$ = {eps2:.4f}')
    plt.tight_layout()
    plt.show()
