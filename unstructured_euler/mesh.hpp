#pragma once
#include <array>
#include <igl/edge_lengths.h>
#include <igl/readMESH.h>
#include <igl/triangle_triangle_adjacency.h>
#include <string>

#include "writer.hpp"

//! The Mesh class is a simple representation of a triangular, unstructered mesh
class Mesh {
  public:
    enum BoundaryType { INTERIOR_EDGE, WING_EDGE, OUTFLOW_EDGE };

    //! Constructs the mesh by reading it from file
    //!
    //! @param filename the name of the file from which to read
    //!
    //! @note assumes the file is MESH formated.
    explicit Mesh(const std::string &filename,
                  const std::function<BoundaryType(const Eigen::Vector2d &,
                                                   const Eigen::Vector2d &)>
                      &classify_boundary) {
        Eigen::MatrixXi tetrahedra;
        igl::readMESH(filename, vertices, tetrahedra, triangles);
        igl::triangle_triangle_adjacency(
            triangles, neighbouringTriangles, neighbouringEdges);

        createAreas();
        createNormals();
        createCellCenters();
        createEdgeCenters();
        classifyBoundaryEdges(classify_boundary);
    }

    //! Gets the number of triangles
    int getNumberOfTriangles() const { return triangles.rows(); }

    //! Gets the indices of the triangle
    Eigen::Vector3i getTriangle(int triangle) const {
        return triangles.row(triangle);
    }

    //! Gets the index of the neighbouring triangle at the given edge,
    //! returns a negative number if there is no neighbour
    int getNeighbour(int triangle, int edge) const {
        return neighbouringTriangles(triangle, edge);
    }

    //! Is neighbour 'triangle' a valid neighbour of 'edge'?
    bool isValidNeighbour(int triangle, int edge) const {
        return getNeighbour(triangle, edge) >= 0;
    }

    //! Gets the normal vector at the given triangle and edge. The normal
    //! vector will point out of the triangle and have length equal to the
    //! length of the edge.
    Eigen::Vector2d getUnitNormal(int triangle, int edge) const {
        return normals.row(triangle * 3 + edge);
    }

    //! Index of the k-th vertex of cell i.
    int getVertexIndex(int i, int k) const { return triangles(i, k); }

    //! Coordinates of the k-th vertex of cell 'i'.
    Eigen::Vector2d getVertex(int i, int k) const {
        int iv = getVertexIndex(i, k);

        Eigen::Vector2d v;
        v << vertices(iv, 0), vertices(iv, 1);
        return v;
    }

    //! Coordinates of the center of the k-th edge of cell 'i'.
    Eigen::Vector2d getEdgeCenter(int i, int k) const {
        return edge_centers.row(i * 3 + k);
    }

    //! Length of the k-th edge of cell i.
    double getEdgeLength(int i, int k) const { return edge_lengths(i, k); }

    //! Coordinates of the center of cell 'i'.
    Eigen::Vector2d getCellCenter(int i) const { return cell_centers.row(i); }

    //! Gets the smallest incircle in the mesh.
    double getMinimumInradius() const {
        double min_inradius = std::numeric_limits<double>::max();

        int n_cells = getNumberOfTriangles();
        for (int i = 0; i < n_cells; ++i) {
            double a = getEdgeLength(i, 0);
            double b = getEdgeLength(i, 1);
            double c = getEdgeLength(i, 2);
            double s = 0.5 * (a + b + c);

            // This is Heron's formula to compute the inradius
            double r = std::sqrt(s * (s - a) * (s - b) * (s - c)) / s;
            min_inradius = std::min(min_inradius, r);
        }

        return min_inradius;
    }

    //! Gets the area of the given triangle
    double getTriangleArea(int triangle) const { return areas(triangle); }

    //! Gets the boundary type of the given edge
    //!
    //! @returns INTERIOR_EDGE if it is not a boundary edge
    //!          WING_EDGE if it is an edge on the wing surface
    //!          OUTFLOW_EDGE if it is an edge on the outflow boundary
    BoundaryType getBoundaryType(int triangle, int edge) const {
        return boundaryTypes[triangle * 3 + edge];
    }

    void writeToFile(const std::string &filename) const {
        writeMatrixToFile(filename + "_vertices.txt", vertices);
        writeMatrixToFile(filename + "_triangles.txt", triangles);
        writeMatrixToFile(filename + "_normals.txt", normals);
        writeMatrixToFile(filename + "_volumes.txt", areas);
        writeMatrixToFile(filename + "_cell_centers.txt", cell_centers);
        writeMatrixToFile(filename + "_edge_centers.txt", edge_centers);
    }

  private:
    void createAreas() {
        int n_cells = triangles.rows();
        areas.resize(n_cells);

        for (int i = 0; i < n_cells; ++i) {
            Eigen::Vector2d a = getVertex(i, 0);
            Eigen::Vector2d b = getVertex(i, 1);
            Eigen::Vector2d c = getVertex(i, 2);

            auto v = b - a;
            auto w = c - a;
            areas(i) = 0.5 * std::abs(v[0] * w[1] - w[0] * v[1]);
        }
    }

    void createNormals() {
        int n_cells = triangles.rows();
        normals.resize(n_cells * 3, 2);
        edge_lengths.resize(n_cells, 3);

        for (int i = 0; i < n_cells; ++i) {
            for (int k = 0; k < 3; ++k) {
                Eigen::Vector2d a = getVertex(i, k % 3);
                Eigen::Vector2d b = getVertex(i, (k + 1) % 3);
                Eigen::Vector2d tau = b - a;

                Eigen::Vector2d normal;
                normal << tau[1], -tau[0];

                edge_lengths(i, k) = normal.norm();
                normals.row(i * 3 + k) = normal / edge_lengths(i, k);
            }
        }
    }

    void createCellCenters() {
        int n_cells = getNumberOfTriangles();
        cell_centers.resize(n_cells, 2);

        for (int i = 0; i < n_cells; ++i) {
            cell_centers.row(i)
                = (getVertex(i, 0) + getVertex(i, 1) + getVertex(i, 2)) / 3.0;
        }
    }

    void createEdgeCenters() {
        int n_cells = getNumberOfTriangles();
        edge_centers.resize(3 * n_cells, 2);

        for (int i = 0; i < n_cells; ++i) {
            for (int k = 0; k < 3; ++k) {
                edge_centers.row(i * 3 + k)
                    = 0.5 * (getVertex(i, k) + getVertex(i, (k + 1) % 3));
            }
        }
    }

    void classifyBoundaryEdges(
        const std::function<BoundaryType(const Eigen::Vector2d &,
                                         const Eigen::Vector2d &)>
            &classify_boundary) {
        boundaryTypes.resize(triangles.rows() * 3);
        for (int triangle = 0; triangle < triangles.rows(); ++triangle) {
            for (int face = 0; face < 3; ++face) {
                int outputIndex = 3 * triangle + face;
                int neighbour = getNeighbour(triangle, face);

                if (neighbour < 0) {
                    auto a = getVertex(triangle, face % 3);
                    auto b = getVertex(triangle, (face + 1) % 3);

                    boundaryTypes[outputIndex] = classify_boundary(a, b);

                } else {
                    boundaryTypes[outputIndex] = INTERIOR_EDGE;
                }
            }
        }
    }

    //! Contains all the vertices of the mesh
    Eigen::MatrixXd vertices;

    //! Contains the triangles of the mesh
    Eigen::MatrixXi triangles;

    //! Contains the edges of the mesh
    Eigen::MatrixXi edges;

    //! Contains the cell-centers.
    Eigen::MatrixXd cell_centers;

    //! For each edge (i, j) it stores the edge_center.
    Eigen::MatrixXd edge_centers;

    //! For each triangle i, and each edge [0,1,2], contains the
    //! index of the neighbouring triangle.
    //!
    //! \note the edges are ordered {0,1}, {1,2} and {2,3}
    Eigen::MatrixXi neighbouringTriangles;

    //! For each triangle i, and each edge [0,1,2], contains the
    //! index of the edge for the neighbouring triangle.
    //!
    //! \note the edges are ordered {0,1}, {1,2} and {2,3}
    Eigen::MatrixXi neighbouringEdges;

    //! For each triangle, for each edge, contains the normal accross the
    //! edge.
    Eigen::MatrixXd normals;

    //! Length of each edge (i, k).
    Eigen::MatrixXd edge_lengths;

    //! Area of each triangle.
    Eigen::VectorXd areas;

    //! For each triangle, for each edge, store the type of the edge
    std::vector<BoundaryType> boundaryTypes;
};

Mesh::BoundaryType classify_naca_boundary(const Eigen::Vector2d &a,
                                          const Eigen::Vector2d &b) {
    if (a.norm() > 10 || b.norm() > 10) {
        return Mesh::BoundaryType::OUTFLOW_EDGE;
    } else {
        return Mesh::BoundaryType::WING_EDGE;
    }
};

Mesh::BoundaryType classify_outflow_boundary(const Eigen::Vector2d &,
                                             const Eigen::Vector2d &) {
    return Mesh::BoundaryType::OUTFLOW_EDGE;
}

Mesh load_naca_mesh(const std::string &filename) {
    return Mesh(filename, classify_naca_boundary);
}

Mesh load_square_mesh(const std::string &filename) {
    return Mesh(filename, classify_outflow_boundary);
}
