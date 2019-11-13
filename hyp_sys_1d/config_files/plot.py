import sys
import json
import glob
import time

import numpy as np
import matplotlib.pyplot as plt

animate_live = False

def load_json(filename):
    with open(filename, "r") as f:
        return json.load(f)


def load_snapshot(filename):
    j = load_json(filename)
    return j["time"], np.array(j["cell_centers"]), np.array(j["density"]), np.array(j["velocity"]), np.array(j["pressure"])


def plot_snapshot(t, x, rho, u, p):
    plt.plot(x, rho, marker="s", linestyle="-", label="$\\rho$")
    plt.plot(x, u, marker="s", linestyle="-", label="$u$")
    plt.plot(x, p, marker="s", linestyle="-", label="$p$")
    plt.title(f"t = {t:0.3f}")


if __name__ == "__main__":
    stem = sys.argv[1]

    files = sorted(glob.glob(f"{stem}*.json"))

    if animate_live:
        plt.ion()
        plt.show()
    else:
        imgs = []

    for fn in files:
        t, x, rho, u, p = load_snapshot(fn)
        plot_snapshot(t, x, rho, u, p)
        plt.legend()
        if animate_live:
            plt.pause(0.05)
        else:
            print(f"Saving {fn}...")
            fn_nofmt = fn[:-5]
            plt.savefig(f"{fn_nofmt}", bbox_inches="tight")
            imgs.append(fn_nofmt)
        plt.clf()


    if not animate_live:
        try:
            import imageio
            read_imgs = [imageio.imread(filename + ".png") for filename in imgs]
            imageio.mimsave(f"{stem}animation.gif", read_imgs, loop=1)
        except ModuleNotFoundError:
            print("Module imageio not found. You might want to try: \n \
             pip install imageio")
