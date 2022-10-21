#include "grid/regular_grid_boundary.hpp"

RegularGridBoundary::RegularGridBoundary(DirectionCode direction, const std::vector<int>& point_indices):
		direction(direction), _point_indices(point_indices){
	std::sort(_point_indices.begin(), _point_indices.end());
}


void RegularGridBoundary::remove_points(const std::vector<int>& indices_to_remove){
	for (int ip: indices_to_remove){
		auto fnd = std::lower_bound(_point_indices.begin(), _point_indices.end(), ip);
		if (fnd != _point_indices.end() && *fnd == ip){
			_point_indices.erase(fnd);
		}
	}
}
