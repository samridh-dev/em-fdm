import h5py
import numpy as np

import os
from multiprocessing import Pool

import matplotlib.pyplot as plt

cmap = 'seismic'
def plot_single_timestep(args):
    t_index, E, B, J, out_dir, n_y = args

    J0 = J[t_index, :, 0::2]
    B_ = B[t_index, :, :]

    E0 = E[t_index, :, 0::2]
    E1 = E[t_index, :, 1::2]

    # Create subplots
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns

    # Plot B_
    axs[0].imshow(B_, aspect='auto', cmap=cmap, origin='lower')
    axs[0].set_title('B')
    axs[0].set_xlabel('X-axis')
    axs[0].set_ylabel('Y-axis')

    # Plot E0
    axs[1].imshow(E0, aspect='auto', cmap=cmap, origin='lower')
    axs[1].set_title('E0')
    axs[1].set_xlabel('X-axis')
    axs[1].set_ylabel('Y-axis')

    # Plot E1
    axs[2].imshow(E1, aspect='auto', cmap=cmap, origin='lower')
    axs[2].set_title('E1')
    axs[2].set_xlabel('X-axis')
    axs[2].set_ylabel('Y-axis')

    plt.tight_layout()
    fout = f"{out_dir}/t{t_index}.png"
    plt.savefig(fout, bbox_inches="tight")
    plt.close(fig)

    print(f"Generated timestep: {t_index}", end="\r")


def create_plots_parallel(h5_file, out_dir="plots"):
    """
    Read HDF5 file and create side-by-side heatmaps in parallel.
    """
    with h5py.File(h5_file, "r") as f:
        E = f["E"][...]  
        B = f["B"][...]
        J = f["J"][...]

    n_t = E.shape[0]
    ny2 = E.shape[2]
    n_y = ny2 // 2

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    args = [(t, E, B, J, out_dir, n_y) for t in range(n_t)]
    with Pool(processes=os.cpu_count()) as pool:
        pool.map(plot_single_timestep, args)

if __name__ == "__main__":
    fname = "out/dat.h5"
    outdir = "out/"
    create_plots_parallel(fname, outdir)
