#include <sstream>
#include <fstream>
#include <numeric>
#include "linear_fem_approximator.hpp"
#include "fe/quadratic_segment.hpp"
#include "fe/linear_segment.hpp"
#include "fe/linear_triangle.hpp"
#include "fe/linear_tetrahedron.hpp"
#include "fe/bilinear_quadrangle.hpp"
#include "debug.hpp"

std::shared_ptr<LinearFemApproximator> LinearFemApproximator::build(std::shared_ptr<AGrid> grid){
	std::shared_ptr<LinearFemApproximator> ret(new LinearFemApproximator(grid));
	ret->initialize();
	return ret;
}

LinearFemApproximator::LinearFemApproximator(std::shared_ptr<AGrid> grid): ANodalFemApproximator(grid){
}

int LinearFemApproximator::n_bases() const {
	return _grid->n_points();
}

Point LinearFemApproximator::node(int inode) const{
	if (inode < _grid->n_points()){
		return _grid->point(inode);
	}
	_THROW_UNREACHABLE_;
}


AElement* LinearFemApproximator::build_element(int icell){
	std::vector<int> cp = _grid->tab_cell_point(icell);
	std::vector<Point> points;

	for (int ipoint: cp){
		points.push_back(_grid->point(ipoint));
	}

	if (_grid->dim==1){
		if (cp.size() == 2){
			return new LinearSegmentElement(cp, points);
		}
	} else if (_grid->dim==2){
		if (cp.size() == 3){
			return new LinearTriangleElement(cp, points);
		} else if (cp.size() == 4){
			return new BilinearQuadrangleElement(cp, points);
		}
	} else if (_grid->dim==3){
		if (cp.size() == 4){
			return new LinearTetrahedronElement(cp, points);
		}
	}

	std::stringstream oss;
	oss << "Failed to find appropriate linear element ";
	oss << "for dim=" << _grid->dim;
	oss << ", npoints=" << cp.size();

	throw std::runtime_error(oss.str());
}

AElement* LinearFemApproximator::build_boundary_element(int iface){
	std::vector<int> cp = _grid->tab_face_point(iface);
	std::vector<Point> points;

	for (int ipoint: cp){
		points.push_back(_grid->point(ipoint));
	}

	std::array<int, 2> cells = _grid->tab_face_cell(iface);
	int icell = (cells[0] >= 0) ? cells[0] : cells[1];
	Point cell_center = _grid->cell_center(icell);
	std::array<double, 4> face_eq = _grid->face_plane(iface);
	double sgn_cc = face_eq[0]*cell_center.x + face_eq[1]*cell_center.y + face_eq[2]*cell_center.z + face_eq[3];
	int direction = sgn_cc > 0 ? -1 : 1;

	if (_grid->dim == 1){
		if (cp.size() == 1){
			return new LinearSegmentBoundaryElement(cp, points);
		}
	} else if (_grid->dim == 2){
		if (cp.size() == 2){
			return new LinearTriangleBoundaryElement(cp, points, direction);
		}
	} else if (_grid->dim == 3){
		if (cp.size() == 3){
			return new LinearTetrahedronBoundaryElement(cp, points, direction);
		}
	}

	std::stringstream oss;
	oss << "Failed to find appropriate linear element ";
	oss << "for dim=" << _grid->dim;
	oss << ", npoints=" << cp.size();

	throw std::runtime_error(oss.str());
}
