#! /usr/bin/env python3

from numpy import *
from pylab import *

import matplotlib.tri
from math import atan2

import sys

vertices = loadtxt(sys.argv[1] + "_vertices.txt")
indices = loadtxt(sys.argv[1] + "_triangles.txt")
normals = loadtxt(sys.argv[1] + "_normals.txt")
edge_centers = loadtxt(sys.argv[1] + "_edge_centers.txt")
cell_centers = loadtxt(sys.argv[1] + "_cell_centers.txt")

values = loadtxt(sys.argv[1] + "_u.txt")

grid = matplotlib.tri.Triangulation(vertices[:, 0], vertices[:, 1], indices)


# # Plot wireframe of the grid.
# triplot(vertices[:, 0], vertices[:, 1], indices)
# show()


mx = values[:, 1]
my = values[:, 2]
rho = values[:, 0]

tripcolor(grid, rho)
colorbar()
# quiver(x,y,nx,ny)
title("rho")
show()


show()
quiver(cell_centers[:, 0], cell_centers[:, 1], mx, my)
title("m")
show()
quiver(cell_centers[:, 0], cell_centers[:, 1], mx / rho, my / rho)
title("u")


tripcolor(grid, values[:, 1] / values[:, 0])
title("mx")
colorbar()
show()

tripcolor(grid, values[:, 2] / values[:, 0])
title("my")
colorbar()
show()


tripcolor(grid, values[:, 3])
title("E")
colorbar()
show()
