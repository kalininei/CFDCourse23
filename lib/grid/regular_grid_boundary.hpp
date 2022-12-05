#ifndef GRID_BOUNDARY_HPP
#define GRID_BOUNDARY_HPP

#include "common.hpp"
#include "grid/agrid_boundary.hpp"

enum struct DirectionCode{
	X_PLUS,
	X_MINUS,
	Y_PLUS,
	Y_MINUS,
	Z_PLUS,
	Z_MINUS
};

class ARegularGrid;

class RegularGridBoundary: public AGridBoundary{
public:
	RegularGridBoundary(DirectionCode direction, const ARegularGrid* grid);
	RegularGridBoundary(DirectionCode direction, std::function<bool(Point)> filter, const ARegularGrid* grid);

	// direction
	const DirectionCode direction;
};


#endif
