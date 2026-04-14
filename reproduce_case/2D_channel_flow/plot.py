import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob

# =========================
# 1. 读取数据文件
# =========================
files = sorted(glob.glob("lbm_output/ux_*.npy"))

# 读第一帧确定尺寸
sample = np.load(files[0])
ny, nx = sample.shape

umax_display = 0.004

# =========================
# 2. 初始化画布
# =========================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

im = ax1.imshow(
    sample,
    origin='lower',
    aspect='auto',
    vmin=0.0,
    vmax=umax_display,
    cmap='jet'
)
plt.colorbar(im, ax=ax1)
ax1.set_title("ux field")

y = np.arange(ny)
x_mid = nx // 2

line_lbm, = ax2.plot(sample[:, x_mid], y, 'o-', label='LBM')
line_theory, = ax2.plot(np.zeros_like(y), y, '-', label='Theory')

ax2.set_xlim(0, umax_display)
ax2.set_ylim(0, ny - 1)
ax2.legend()
ax2.grid(True)

time_text = ax2.text(0.05, 0.95, "", transform=ax2.transAxes, va='top')

# =========================
# 3. 更新函数
# =========================
def update(i):
    ux = np.load(files[i])

    im.set_array(ux)

    u_profile = ux[:, x_mid]
    line_lbm.set_data(u_profile, y)

    # 理论抛物线
    y_phys = y / (ny - 1)
    u_theory = y_phys * (1 - y_phys)
    if np.max(u_profile) > 0:
        u_theory *= np.max(u_profile) / np.max(u_theory)

    line_theory.set_data(u_theory, y)

    time_text.set_text(f"frame = {i}")

    return im, line_lbm, line_theory, time_text

# =========================
# 4. 动画
# =========================
ani = FuncAnimation(fig, update, frames=len(files), interval=50)

plt.tight_layout()
plt.show()

# 保存视频（可选）
ani.save("visualization.mp4", writer="ffmpeg", fps=20)