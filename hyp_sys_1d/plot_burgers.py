import json
import glob

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import LogLocator
from matplotlib import rc
from matplotlib import rcParams
font = {'family' : 'Dejavu Sans',
        'weight' : 'normal',
        'size'   : 22}
rc('font', **font)
rcParams['lines.linewidth'] = 4
rcParams['lines.markersize'] = 12
rcParams['markers.fillstyle'] = 'none'

def load_json(filename):
    with open(filename, "r") as f:
        return json.load(f)

def load_snapshot(filename):
    j = load_json(filename)
    return j["time"], np.array(j["cell_centers"]), np.array(j["data"])

def plot_snapshot(t, x, u):
    plt.plot(x, u)
    plt.title(f"t = {t:0.6f}")
    plt.show()

if __name__ == "__main__":

    files = sorted(glob.glob("output/burgers_fvm_rarefaction*.json"))

    for fn in files:
        t, x, u = load_snapshot(fn)
        plot_snapshot(t, x, u) 
