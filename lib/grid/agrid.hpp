#ifndef GRID_AGRID_HPP
#define GRID_AGRID_HPP

#include "common.hpp"
#include "geom.hpp"
#include "grid/agrid_boundary.hpp"

class AGrid{
public:
	virtual ~AGrid() = default;
	AGrid(int ndim);

	// geometric dimension
	const int dim;

	// sizes
	virtual int n_points() const = 0;
	virtual int n_faces() const = 0;

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

	// face->point connectivity
	// points are ordered in counter clockwise direction
	virtual std::vector<int> tab_face_point(int iface) const = 0;
	
	// cell->faces connectivity
	virtual std::vector<int> tab_cell_face(int icell) const = 0;

	// face->cell connectivity
	// [0] - left cell: 1d - cell with lower x
	//                  2d - cell to the left hand side of the directed face
	//                  3d - top cell if we look from the side where face points have cc ordering
	// [1] - right cell
	//
	// If there is no left or right cell, returns -1 on respective places
	virtual std::array<int, 2> tab_face_cell(int iface) const = 0;

	// cell->point connectivity
	virtual std::vector<int> tab_cell_point(int icell) const = 0;

	// ==== boundaries
	// passed boundary types
	std::vector<int> btypes() const;

	// returns boundary for the specified type
	const AGridBoundary& boundary(int btype) const;

	// define boundary
	virtual void define_boundary(int btype, std::shared_ptr<AGridBoundary> boundary);
private:
	std::map<int, std::shared_ptr<AGridBoundary>> _boundary;
};
#endif
