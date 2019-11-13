#! /usr/bin/env python3

import scipy.interpolate
import numpy as np
import sys
import json

GAMMA = 1.4


class Grid:
    def __init__(self, vertices, triangles, cell_centers, volumes):
        self.vertices = vertices
        self.triangles = triangles
        self.cell_centers = cell_centers
        self.volumes = volumes


def load_grid(basename):
    vertices = np.loadtxt(basename + "_vertices.txt")
    triangles = np.loadtxt(basename + "_triangles.txt")
    cell_centers = np.loadtxt(basename + "_cell_centers.txt")
    volumes = np.loadtxt(basename + "_volumes.txt")

    return Grid(vertices, triangles, cell_centers, volumes)


def load_soln(basename):
    return np.loadtxt(basename + "_u.txt")


class ReferenceSoln:
    def __init__(self, grid, u_ref):
        """Args:
             grid: cell_centers. shape = (n_cells, 2)
             u_ref: values of the reference solution. shape = (n_cells, ...)
        """

        x_ref = grid.cell_centers

        self.linear_interpolate = scipy.interpolate.LinearNDInterpolator(
            x_ref, u_ref, fill_value=np.nan
        )

        self.nearest_interpolate = scipy.interpolate.NearestNDInterpolator(x_ref, u_ref)

    def __call__(self, x):
        """Returns the reference solution at x.

        Args:
            x: shape = (n_coarse_cells, 2)
        """

        u_ref = self.linear_interpolate(x)

        # If we're unlucky some points don't lie in the convex hull of the
        # cell-centers of the reference grid.
        # In this case we would like to fall-back to nearest neighbour interpolation.

        I = np.all(np.isnan(u_ref), axis=1)
        u_ref[I, :] = self.nearest_interpolate(x[I, :])

        return u_ref

    def l1_error(self, grid_coarse, u_coarse):
        """Compute the L1-error against the reference solution.

        Args:
            grid_coarse: the coarse grid.
            u_coarse: approximation of 'u' on the coarse grid.
        """

        vol = np.reshape(grid_coarse.volumes, (-1, 1))
        x_coarse = grid_coarse.cell_centers

        return np.sum(vol * np.abs(self(x_coarse) - u_coarse), axis=0)


def load_reference_soln(basename):
    grid_ref = load_grid(basename)
    u_ref = load_soln(basename)

    return ReferenceSoln(grid_ref, u_ref)


def compute_vortex_error(ref_soln, basename):
    grid_coarse = load_grid(basename)
    u_coarse = load_soln(basename)

    return ref_soln.l1_error(grid_coarse, u_coarse)


if __name__ == "__main__":
    if len(sys.argv) <= 4:
        print(f"Usage: {sys.argv[0]} APPROX_1 APPROX_2 [APPROX_3, ...] REFERENCE")
        print("")
        print("Where 'APPROX_0' is the basename of the experiment of the")
        print("coarsest level and REFERENCE is the basename on the finest level.")
        print("")
        print("The refinement from one level should be '2', e.g. approx 4x cells.")

    ref_soln = load_reference_soln(sys.argv[-1])

    errors = np.array([compute_vortex_error(ref_soln, f) for f in sys.argv[1:-1]])
    rates = (np.log(errors[:-1]) - np.log(errors[1:])) / np.log(2.0)

    print(errors)
    print(rates)
