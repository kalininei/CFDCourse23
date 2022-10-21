#ifndef GRID_BOUNDARY_HPP
#define GRID_BOUNDARY_HPP

#include "common.hpp"

enum struct DirectionCode{
	X_PLUS,
	X_MINUS,
	Y_PLUS,
	Y_MINUS,
	Z_PLUS,
	Z_MINUS
};

class RegularGridBoundary{
public:
	RegularGridBoundary(DirectionCode direction, const std::vector<int>& point_indices);

	// direction
	const DirectionCode direction;
	const std::vector<int>& point_indices() const { return _point_indices; }

	void remove_points(const std::vector<int>& indices_to_remove);
private:
	// boundary points indices
	std::vector<int> _point_indices;
};


#endif
