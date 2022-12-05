#include "agrid_boundary.hpp"
#include "agrid.hpp"

namespace{

std::vector<int> build_point_indices(const std::vector<int>& face_indices, const AGrid* grid){
	std::vector<int> ret;

	for (int iface: face_indices){
		std::vector<int> ipoints = grid->tab_face_point(iface);
		ret.insert(ret.end(), ipoints.begin(), ipoints.end());
	}

	std::sort(ret.begin(), ret.end());
	ret.resize(std::unique(ret.begin(), ret.end()) - ret.begin());

	return ret;
}

std::vector<int> build_cell_indices(const std::vector<int>& face_indices, const AGrid* grid){
	std::vector<int> ret;

	for (int iface: face_indices){
		std::array<int, 2> icells = grid->tab_face_cell(iface);
		if (icells[0] >= 0){
			ret.push_back(icells[0]);
		} else {
			ret.push_back(icells[1]);
		}
	}

	return ret;
}

std::vector<int> apply_filter(const std::vector<int>& face_indices, std::function<bool(Point)> filter, const AGrid* grid){
	std::vector<int> ret;

	for (int iface: face_indices){
		if (filter(grid->face_center(iface))){
			ret.push_back(iface);
		}
	}

	if (ret.size() == 0){
		throw std::runtime_error("No boundary faces with the given filter are detected");
	}

	return ret;
}

}

AGridBoundary::AGridBoundary(const std::vector<int>& face_indices, const AGrid* grid)
	: _point_indices(build_point_indices(face_indices, grid)),
	  _face_indices(face_indices),
	  _cell_indices(build_cell_indices(face_indices, grid)){
}

AGridBoundary::AGridBoundary(const std::vector<int>& face_indices, std::function<bool(Point)> filter, const AGrid* grid)
	: AGridBoundary(apply_filter(face_indices, filter, grid), grid){
}

int AGridBoundary::n_points() const{
	return (int)_point_indices.size();
}
int AGridBoundary::n_faces() const{
	return (int)_face_indices.size();
}

const std::vector<int>& AGridBoundary::grid_point_indices() const{
	return _point_indices;
}

const std::vector<int>& AGridBoundary::grid_face_indices() const{
	return _face_indices;
}

const std::vector<int>& AGridBoundary::grid_cell_indices() const{
	return _cell_indices;
}
