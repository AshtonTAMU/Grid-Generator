import numpy as np
import matplotlib.pyplot as plt

# --- USER SETTINGS ---
DAT_FILE = "grid.dat"
FIGSIZE  = (10, 7)
LINECOLOR = "steelblue"
LINEWIDTH = 0.6
TITLE     = "Laplace Smoothed Structured Grid"
# ----------------------

def read_tecplot_grid(filename):
    """Read a Tecplot ASCII POINT-format zone with I,J dimensions."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    ni = nj = None
    data_start = 0

    for idx, line in enumerate(lines):
        low = line.lower()
        if 'zone' in low and 'i=' in low:
            # parse  I=41, J=31
            parts = line.replace(',', ' ').replace('=', ' ').split()
            for k, p in enumerate(parts):
                if p.upper() == 'I':
                    ni = int(parts[k+1])
                if p.upper() == 'J':
                    nj = int(parts[k+1])
            data_start = idx + 1
            break

    if ni is None or nj is None:
        raise ValueError("Could not parse I,J dimensions from ZONE header.")

    raw = []
    for line in lines[data_start:]:
        vals = line.split()
        if len(vals) >= 2:
            raw.append((float(vals[0]), float(vals[1])))

    if len(raw) < ni * nj:
        raise ValueError(f"Expected {ni*nj} points, got {len(raw)}.")

    x = np.array([p[0] for p in raw[:ni*nj]]).reshape(nj, ni)
    y = np.array([p[1] for p in raw[:ni*nj]]).reshape(nj, ni)
    return x, y, ni, nj


def plot_grid(x, y, ni, nj):
    fig, ax = plt.subplots(figsize=FIGSIZE)

    # Draw i-lines (constant j index → horizontal-ish lines)
    for j in range(nj):
        ax.plot(x[j, :], y[j, :], color=LINECOLOR, linewidth=LINEWIDTH)

    # Draw j-lines (constant i index → vertical-ish lines)
    for i in range(ni):
        ax.plot(x[:, i], y[:, i], color=LINECOLOR, linewidth=LINEWIDTH)

    ax.set_aspect('equal')
    ax.set_title(TITLE, fontsize=13)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    plt.tight_layout()
    plt.savefig("grid_plot.png", dpi=150)
    print("Saved grid_plot.png")
    plt.show()


if __name__ == "__main__":
    x, y, ni, nj = read_tecplot_grid(DAT_FILE)
    print(f"Grid loaded: {ni} x {nj} points")
    plot_grid(x, y, ni, nj)
