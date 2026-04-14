import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# =========================
# 1. 用户参数
# =========================

parameter = 'tau1_R25'
mode = "density"             # 可选: "density", "umag", "density_quiver"

folder = f"output_{parameter}"   # 改成你的输出文件夹
interval = 200               # 动画帧间隔(ms)
stride = 10                   # 箭头稀疏步长，仅 density_quiver 模式用
save_gif = True              # 是否保存 gif
gif_name = f"droplet_animation_{mode}_{folder}.gif"

# =========================
# 2. 工具函数
# =========================
def extract_step(filename, prefix):
    m = re.match(rf"{prefix}(\d+)\.dat", filename)
    return int(m.group(1)) if m else None

def load_dat(filepath):
    return np.loadtxt(filepath)

def get_sorted_files(folder, prefix):
    files = []
    for f in os.listdir(folder):
        step = extract_step(f, prefix)
        if step is not None:
            files.append((step, os.path.join(folder, f)))
    files.sort(key=lambda x: x[0])
    return files

# =========================
# 3. 收集文件
# =========================
density_files = get_sorted_files(folder, "density")
velx_files = get_sorted_files(folder, "velx")
vely_files = get_sorted_files(folder, "vely")

if len(density_files) == 0:
    raise FileNotFoundError("没有找到 density*.dat 文件")

if mode in ["umag", "density_quiver"]:
    if len(velx_files) == 0 or len(vely_files) == 0:
        raise FileNotFoundError("需要 velx*.dat 和 vely*.dat 文件")

# 取共同步数，避免文件不齐
if mode == "density":
    steps = [s for s, _ in density_files]
else:
    density_dict = {s: f for s, f in density_files}
    velx_dict = {s: f for s, f in velx_files}
    vely_dict = {s: f for s, f in vely_files}
    common_steps = sorted(set(density_dict) & set(velx_dict) & set(vely_dict))
    if len(common_steps) == 0:
        raise ValueError("density / velx / vely 没有共同时间步")
    density_files = [(s, density_dict[s]) for s in common_steps]
    velx_files = [(s, velx_dict[s]) for s in common_steps]
    vely_files = [(s, vely_dict[s]) for s in common_steps]
    steps = common_steps

# =========================
# 4. 预读取数据
# =========================
densities = [load_dat(f) for _, f in density_files]

if mode in ["umag", "density_quiver"]:
    velxs = [load_dat(f) for _, f in velx_files]
    velys = [load_dat(f) for _, f in vely_files]
    umags = [np.sqrt(vx**2 + vy**2) for vx, vy in zip(velxs, velys)]

ny, nx = densities[0].shape

# 统一颜色范围，避免动画闪烁
rho_min = min(np.min(r) for r in densities)
rho_max = max(np.max(r) for r in densities)

if mode in ["umag", "density_quiver"]:
    u_min = min(np.min(u) for u in umags)
    u_max = max(np.max(u) for u in umags)

# =========================
# 5. 建图
# =========================
fig, ax = plt.subplots(figsize=(6, 6))

if mode == "density":
    im = ax.imshow(densities[0], origin="lower", cmap="jet",
                   vmin=rho_min, vmax=rho_max, animated=True)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Density")
    title = ax.set_title(f"Density, step = {steps[0]}")

    def update(frame):
        im.set_array(densities[frame])
        title.set_text(f"Density, step = {steps[frame]}")
        return [im, title]

elif mode == "umag":
    im = ax.imshow(umags[0], origin="lower", cmap="jet",
                   vmin=u_min, vmax=u_max, animated=True)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("|u|")
    title = ax.set_title(f"Velocity Magnitude, step = {steps[0]}")

    def update(frame):
        im.set_array(umags[frame])
        title.set_text(f"Velocity Magnitude, step = {steps[frame]}")
        return [im, title]

elif mode == "density_quiver":
    im = ax.imshow(densities[0], origin="lower", cmap="jet",
                   vmin=rho_min, vmax=rho_max, animated=True)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Density")

    y = np.arange(0, ny, stride)
    x = np.arange(0, nx, stride)
    X, Y = np.meshgrid(x, y)

    q = ax.quiver(
        X, Y,
        velxs[0][::stride, ::stride],
        velys[0][::stride, ::stride],
        color="white",
        pivot="mid"
    )
    title = ax.set_title(f"Density + Velocity, step = {steps[0]}")

    def update(frame):
        im.set_array(densities[frame])
        q.set_UVC(
            velxs[frame][::stride, ::stride],
            velys[frame][::stride, ::stride]
        )
        title.set_text(f"Density + Velocity, step = {steps[frame]}")
        return [im, q, title]

else:
    raise ValueError("mode 只能是: density, umag, density_quiver")

ax.set_xlabel("x")
ax.set_ylabel("y")

ani = animation.FuncAnimation(
    fig, update, frames=len(steps), interval=interval, blit=False
)

plt.tight_layout()

if save_gif:
    ani.save(gif_name, writer="pillow")
    print(f"动画已保存为: {gif_name}")

plt.show()