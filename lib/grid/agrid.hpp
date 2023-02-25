/// agrid.hpp
#ifndef GRID_AGRID_HPP
#define GRID_AGRID_HPP

#include "common.hpp"
#include "geom.hpp"
#include "grid/agrid_boundary.hpp"

/**
 * @brief Abstract grid type
 *
 * Represents basic grid interface.
 */
class AGrid{
public:
	/// virtual destructor
	virtual ~AGrid() = default;
	AGrid(int ndim);

	// geometric dimension
	const int dim;

	// sizes
	virtual int n_points() const = 0;
	virtual int n_faces() const = 0;
	virtual int n_cells() const = 0;
	virtual int n_boundary_faces() const;
	virtual int n_internal_faces() const;

	// point by the global index
	virtual Point point(int ipoint) const = 0;

	// center of the face
	virtual Point face_center(int iface) const;

	// center of the cell
	virtual Point cell_center(int icell) const;

	// face area
	virtual double face_area(int iface) const = 0;

	// cell volume
	virtual double cell_volume(int icell) const = 0;

	// face normal
	virtual Vector face_normal(int iface) const = 0;

	// list of boundary faces indices
	virtual std::vector<int> boundary_faces() const;

	// list of internal faces indices
	virtual std::vector<int> internal_faces() const;

	// face->point connectivity
	// points are ordered in counter clockwise direction
	virtual std::vector<int> tab_face_point(int iface) const = 0;
	
	// cell->faces connectivity
	virtual std::vector<int> tab_cell_face(int icell) const;

	// cell->cell connectivity
	virtual std::vector<int> tab_cell_cell(int icell) const;

	// face->cell connectivity
	// [0] - left cell: 1d - cell with lower x
	//                  2d - cell to the left hand side of the directed face
	//                  3d - top cell if we look from the side where face points have cc ordering
	// [1] - right cell
	//
	// If there is no left or right cell, returns -1 on respective places
	virtual std::array<int, 2> tab_face_cell(int iface) const = 0;

	// cell->point connectivity
	virtual std::vector<int> tab_cell_point(int icell) const;

	// ==== boundaries
	// passed boundary types
	std::vector<int> btypes() const;

	// returns boundary for the specified type
	const AGridBoundary& boundary(int btype) const;

	// face btypes for all faces
	std::vector<int> face_btypes() const;
	int get_boundary_index_for_face(int iface) const;

	// define boundary
	void define_boundary(int btype, std::function<bool(Point)> face_center_filter);

	// ==== info
	// returns index of the cell which contains given point or -1 if point is out
	// of the grid
	virtual int find_cell_index(Point p) const;

	// returns A, B, C, D coefficients for plane equation in the form: Ax + By + Cz + D = 0
	std::array<double, 4> face_plane(int iface) const;

	// vtk converters
	virtual std::vector<std::vector<int>> vtk_cell_array() const;
	virtual std::vector<int> vtk_cell_types() const;
	virtual std::vector<std::vector<int>> vtk_boundary_face_array() const;
	virtual std::vector<int> vtk_boundary_face_types() const;

protected:
	void define_boundary(int btype, std::shared_ptr<AGridBoundary> boundary);

private:
	struct Cache{
		std::vector<std::vector<int>> tab_cell_point;
		std::vector<std::vector<int>> tab_cell_face;
		std::vector<std::vector<int>> tab_cell_cell;
		std::vector<int> boundary_faces;
		std::vector<int> internal_faces;
		std::map<int, int> grid_to_bnd_face;

		int n_boundary_faces = -1;
		int n_internal_faces = -1;
	};
	mutable Cache _cache;
	std::map<int, std::shared_ptr<AGridBoundary>> _boundary;
};
#endif
