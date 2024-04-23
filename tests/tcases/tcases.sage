import numpy as np
import pandas as pd
#import flin
from datetime import datetime
import matplotlib.pyplot as plt 
import matplotlib
from matplotlib.gridspec import GridSpec
import scipy.special

def fmt(x):
    if isinstance(x, float):
        m, n = f"{x:+.14E}".split("E")
        return f"{float(m):+.14f}E{int(n):+04d}"
    elif isinstance(x, complex):
        return f"{fmt(x.real)} {fmt(x.imag)}"

def compute_r1(alpha, eps):
    assert (alpha > 0 and alpha < 1)
    C0 = (1.3 ** (1 - alpha))/(np.pi * np.sin(np.pi * alpha))
    return (-2.0 * np.log(eps/C0))**alpha

def rnd_circ(condition, r = 15.0, nsamples = 1):
    res, count = [], 0
    while count < nsamples:
        U1, U2 = np.random.uniform(size = 2)
        x = r * np.sqrt(U2) * np.cos(2 * np.pi * U1)
        y = r * np.sqrt(U2) * np.sin(2 * np.pi * U1)
        z = x + 1j * y
        if condition(z):
            res.append(z)
            count += 1
    return np.array(res)

def cross(phi1, phi2):
    return np.sin(phi1 - phi2)
    
def are_clockwise(phi1, phi2, closed = False):
    if closed:
      return cross(phi1, phi2) >= 0
    else:
        return cross(phi1, phi2) > 0

def wedge(z, phi1, phi2, closed = False):
    th = np.angle(z)
    if are_clockwise(phi1, phi2, closed): 
        return not wedge(z, phi2, phi1, closed)
    else: 
        return (not are_clockwise(phi1, th, closed)) and are_clockwise(phi2, th, closed)
        
def create_partition(α = 0.5, ε = 1.0e-15, r0 = 0.95, rmax = 25.0, nsamples = 10000):
    """
    α in (0, 1)
    """
    π = np.pi
    r1 = compute_r1(α, ε)
    δ, τ = π * α / 8.0, min(π * α / 8.0, π * (1.0 + α) / 2.0)
    
    disk = rnd_circ(lambda z: np.abs(z) < rmax, r = 2*r1, nsamples = nsamples)

    in_region_G = [
        lambda z: np.abs(z) <= r0,
        lambda z: np.abs(z) > r1 and wedge(z, -π * α + δ, +π * α - δ),
        lambda z: np.abs(z) > r1 and wedge(z, +π * α + τ, -π * α - τ),
        lambda z: np.abs(z) > r1 and wedge(z, +π * α - δ, +π * α + τ, True),
        lambda z: np.abs(z) > r1 and wedge(z, -π * α - τ, -π * α + δ, True),
        lambda z: np.abs(z) < r1 and np.abs(z) > r0 and wedge(z, -5.0 * π * α / 6.0, +5.0 * π * α / 6.0, True),
        lambda z: np.abs(z) < r1 and np.abs(z) > r0 and wedge(z, +5.0 * π * α / 6.0, -5.0 * π * α / 6.0)
    ]

    tbl = []
    for i in range(7):
        G = np.array(list(filter(in_region_G[i], disk)))
        tbl.append(pd.DataFrame({
            "alpha": α,
            "re_z": G.real,
            "im_z": G.imag,
            "region": i}))
    return pd.concat(tbl)
    
def plot_partition(ax, α = 0.5, ε = 1.0e-15, r0 = 0.95, rmax = 25.0, size = 0.95):
    π = np.pi
    df = create_partition(α, ε, r0, rmax)
    ax.set_title(f"$\\alpha = {α}$")
    ax.set_aspect('equal')
    with plt.style.context({"axes.prop_cycle" : plt.cycler("color", plt.cm.Dark2.colors)}):
        for i in range(7):
            ax.scatter(df[df.region == i].re_z, df[df.region == i].im_z, s = size, label = f"$G_{i}$")
        ax.axline((0, 0), slope = π * α, color = "#414141", linestyle='--')
        ax.axline((0, 0), slope = -π * α, color = "#414141", linestyle='--')
        # Plot legend.
        lgnd = ax.legend(bbox_to_anchor=(1.0, 1.01))
        for i in lgnd.legend_handles:
            i._sizes = [40]
        return ax

print("NumPy version:", np.__version__)
print("Pandas version:", pd.__version__)
print("Matplotlib version:", matplotlib.__version__)
print("Last Updated:", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

