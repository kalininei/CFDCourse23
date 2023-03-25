#include <sstream>
#include <fstream>
#include <numeric>
#include "quadratic_fem_approximator.hpp"
#include "fe/quadratic_segment.hpp"
#include "debug.hpp"

std::shared_ptr<QuadraticFemApproximator> QuadraticFemApproximator::build(std::shared_ptr<AGrid> grid){
	std::shared_ptr<QuadraticFemApproximator> ret;
	if (grid->dim == 1){
		ret.reset(new QuadraticFemApproximator_1D(grid));
	} else {
		_THROW_NOT_IMP_;
	}
	ret->initialize();
	return ret;
}

QuadraticFemApproximator::QuadraticFemApproximator(std::shared_ptr<AGrid> grid): ANodalFemApproximator(grid){
}

// ========================= 1D
QuadraticFemApproximator_1D::QuadraticFemApproximator_1D(std::shared_ptr<AGrid> grid): QuadraticFemApproximator(grid){}

int QuadraticFemApproximator_1D::n_bases() const {
	return _grid->n_points() + _grid->n_cells();
}

Point QuadraticFemApproximator_1D::node(int inode) const{
	if (inode < _grid->n_points()){
		return _grid->point(inode);
	} else {
		int icell = inode - _grid->n_points();
		if (icell < _grid->n_cells()){
			return _grid->cell_center(icell);
		}
	}
	_THROW_UNREACHABLE_;
}

AElement* QuadraticFemApproximator_1D::build_element(int icell){
	std::vector<int> cp = _grid->tab_cell_point(icell);
	std::vector<Point> points;

	for (int ipoint: cp){
		points.push_back(_grid->point(ipoint));
	}

	cp.push_back(_grid->n_points() + icell);
	return new QuadraticSegmentElement(cp, points);
}

AElement* QuadraticFemApproximator_1D::build_boundary_element(int iface){
	std::vector<int> cp = _grid->tab_face_point(iface);
	std::vector<Point> points;

	for (int ipoint: cp){
		points.push_back(_grid->point(ipoint));
	}

	return new QuadraticSegmentBoundaryElement(cp, points);
}
