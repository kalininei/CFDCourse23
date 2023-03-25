#include <sstream>
#include <fstream>
#include <numeric>
#include "cubic_fem_approximator.hpp"
#include "fe/cubic_segment.hpp"
#include "debug.hpp"

std::shared_ptr<CubicFemApproximator> CubicFemApproximator::build(std::shared_ptr<AGrid> grid){
	std::shared_ptr<CubicFemApproximator> ret;
	if (grid->dim == 1){
		ret.reset(new CubicFemApproximator_1D(grid));
	} else {
		_THROW_NOT_IMP_;
	}
	ret->initialize();
	return ret;
}

CubicFemApproximator::CubicFemApproximator(std::shared_ptr<AGrid> grid): ANodalFemApproximator(grid){
}

// ========================= 1D
CubicFemApproximator_1D::CubicFemApproximator_1D(std::shared_ptr<AGrid> grid): CubicFemApproximator(grid){}

int CubicFemApproximator_1D::n_bases() const {
	return _grid->n_points() + 2*_grid->n_cells();
}

Point CubicFemApproximator_1D::node(int inode) const{
	if (inode < _grid->n_points()){
		return _grid->point(inode);
	} else {
		int i2 = inode - _grid->n_points();
		if (i2 < 2*_grid->n_cells()){
			int icell = i2 / 2;
			int n = i2 % 2 + 1;
			double len3 = _grid->cell_volume(icell)/3;
			return _grid->point(icell) + Point(len3*n);
		}
	}
	_THROW_UNREACHABLE_;
}

AElement* CubicFemApproximator_1D::build_element(int icell){
	std::vector<int> cp = _grid->tab_cell_point(icell);
	std::vector<Point> points;

	for (int ipoint: cp){
		points.push_back(_grid->point(ipoint));
	}

	cp.push_back(_grid->n_points() + 2*icell);
	cp.push_back(_grid->n_points() + 2*icell+1);
	return new CubicSegmentElement(cp, points);
}

AElement* CubicFemApproximator_1D::build_boundary_element(int iface){
	std::vector<int> cp = _grid->tab_face_point(iface);
	std::vector<Point> points;

	for (int ipoint: cp){
		points.push_back(_grid->point(ipoint));
	}

	return new CubicSegmentBoundaryElement(cp, points);
}
