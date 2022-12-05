#ifndef GRID_AGRID_BOUNDARY_HPP
#define GRID_AGRID_BOUNDARY_HPP

#include "common.hpp"
#include "geom.hpp"

class AGrid;

class AGridBoundary{
public:
	virtual ~AGridBoundary() = default;
	AGridBoundary(const std::vector<int>& face_indices, const AGrid* grid);
	AGridBoundary(const std::vector<int>& face_indices, std::function<bool(Point)> filter, const AGrid* grid);

	int n_points() const;
	int n_faces() const;

	// List of boundary points for this boundary
	const std::vector<int>& grid_point_indices() const;

	// List of boundary faces for this boundary
	const std::vector<int>& grid_face_indices() const;

	// List of cells adjacent to the boundary faces.
	// Returned vector corresponds to faces indices from AGridBoundary::grid_face_indices()
	const std::vector<int>& grid_cell_indices() const;
private:
	const AGrid* grid;

	const std::vector<int> _point_indices;
	const std::vector<int> _face_indices;
	const std::vector<int> _cell_indices;
};




#endif
