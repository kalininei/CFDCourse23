#include <sstream>
#include <fstream>
#include <numeric>
#include "linear_fem_approximator.hpp"
#include "fe/linear_segment.hpp"
#include "fe/linear_triangle.hpp"
#include "fe/linear_tetrahedron.hpp"

std::shared_ptr<LinearFemApproximator> LinearFemApproximator::build(std::shared_ptr<AGrid> grid){
	return std::shared_ptr<LinearFemApproximator>(new LinearFemApproximator(grid));
}

LinearFemApproximator::LinearFemApproximator(std::shared_ptr<AGrid> grid): _grid(grid){
	for (int icell=0; icell < _grid->n_cells(); ++icell){
		_elements.emplace_back(build_element(icell));
	}

	for (int iface: _grid->boundary_faces()){
		_boundary_elements[iface] = std::shared_ptr<AElement>(build_boundary_element(iface));
	}
}

int LinearFemApproximator::n_bases() const {
	return _grid->n_points();
}

int LinearFemApproximator::n_elements() const {
	return (int)_elements.size();
}

std::vector<double> LinearFemApproximator::approximate(std::function<double(Point)> func) const{
	std::vector<double> ret(n_bases());

	for (int i=0; i<_grid->n_points(); ++i){
		ret[i] = func(_grid->point(i));
	}

	return ret;
}

CsrStencil LinearFemApproximator::_build_stencil() const{
	std::vector<std::set<int>> s(n_bases());

	for (int i=0; i<_grid->n_points(); ++i){
		std::vector<int> pp = _grid->tab_point_point(i);
		s[i].insert(i);
		s[i].insert(pp.begin(), pp.end());
	}

	return CsrStencil::build(s);
}

std::vector<double> LinearFemApproximator::mass() const{
	std::vector<double> ret(stencil().n_nonzero(), 0);

	for (int i=0; i<n_elements(); ++i){
		std::vector<double> local_mass = _elements[i]->mass();
		add_local_matrix(1, local_mass, _elements[i].get(), ret);
	}

	return ret;
}

std::vector<double> LinearFemApproximator::stiff() const{
	std::vector<double> ret(stencil().n_nonzero(), 0);

	for (int i=0; i<n_elements(); ++i){
		std::vector<double> local_stiff = _elements[i]->stiff();
		add_local_matrix(1, local_stiff, _elements[i].get(), ret);
	}

	return ret;
}

std::vector<double> LinearFemApproximator::_build_load_vector() const{
	std::vector<double> ret(stencil().n_rows(), 0);

	for (int i=0; i<n_elements(); ++i){
		std::vector<double> local_load = _elements[i]->load();
		add_local_vector(1, local_load, _elements[i].get(), ret);
	}

	return ret;
}

void LinearFemApproximator::apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const{
	const AGridBoundary& bnd = _grid->boundary(ibnd);
	for (int iface: bnd.grid_face_indices()){
		const AElement* el = _boundary_elements.find(iface)->second.get();
		std::vector<double> locmass = el->mass();
		std::vector<double> locvec(el->n_bases(), 0);

		for (int i=0; i<el->n_bases(); ++i){
			for (int j=0; j<el->n_bases(); ++j){
				int icol = el->i_basis(j);
				locvec[i] += locmass[i * el->n_bases() + j] * q_func(_grid->point(icol));
			}
		}

		add_local_vector(1, locvec, el, rhs);
	}
}


void LinearFemApproximator::add_local_matrix(double coef, const std::vector<double>& lmat, const AElement* elem, std::vector<double>& gmat) const{
	int k = 0;
	for (int j=0; j<elem->n_bases(); ++j){
		int irow = elem->i_basis(j);
		for (int i=0; i<elem->n_bases(); ++i){
			int icol = elem->i_basis(i);
			int addr = stencil().addr_index(irow, icol);
			gmat[addr] += coef*lmat[k++];
		}
	}
}

void LinearFemApproximator::add_local_vector(double coef, const std::vector<double>& lvec, const AElement* elem, std::vector<double>& gvec) const{
	int k = 0;
	for (int j=0; j<elem->n_bases(); ++j){
		int irow = elem->i_basis(j);
		gvec[irow] += coef*lvec[k++];
	}
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

std::map<int, std::vector<std::pair<int, Point>>> LinearFemApproximator::_build_boundary_bases() const{
	std::map<int, std::vector<std::pair<int, Point>>> ret;

	// sort btypes in reverse order to give higher priority to large btype values
	std::vector<int> btypes = _grid->btypes();
	std::sort(btypes.begin(), btypes.end());
	std::reverse(btypes.begin(), btypes.end());

	// using std::set to prevent duplicating of points in different bc faces
	std::set<int> used_points;

	for (int btype: btypes){
		std::vector<std::pair<int, Point>>& vec = ret[btype];
		const AGridBoundary& bnd = _grid->boundary(btype);
		for (int ipoint: bnd.grid_point_indices()){
			std::pair<std::set<int>::const_iterator, bool> ires = used_points.insert(ipoint);
			// if point was not used before
			if (ires.second){
				vec.push_back({ipoint, _grid->point(ipoint)});
			}
		}
	}

	return ret;
}

void LinearFemApproximator::vtk_save_scalars(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const{
	std::ofstream ofs(filepath);

	// header
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "fvm output" << std::endl;
	ofs << "ASCII" << std::endl;

	// grid
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << "POINTS " << _grid->n_points() << " double" << std::endl;
	for (int ipoint=0; ipoint<_grid->n_points(); ++ipoint){
		auto point = _grid->point(ipoint);
		ofs << point.x << " " << point.y << " " << point.z << std::endl;
	}

	int nc = _grid->n_cells() + _grid->n_boundary_faces();

	std::vector<std::vector<int>> cell_array = _grid->vtk_cell_array();
	std::vector<int> cell_types = _grid->vtk_cell_types();
	std::vector<std::vector<int>> bnd_cell_array = _grid->vtk_boundary_face_array();
	std::vector<int> bnd_cell_types = _grid->vtk_boundary_face_types();
	cell_array.insert(cell_array.end(), bnd_cell_array.begin(), bnd_cell_array.end());
	cell_types.insert(cell_types.end(), bnd_cell_types.begin(), bnd_cell_types.end());

	int nfull = std::accumulate(
		cell_array.begin(), cell_array.end(), 0,
		[](int ret, const std::vector<int>& v)->int{ return ret + int(v.size()) + 1; });
	ofs << "CELLS " << nc << " " << nfull << std::endl;
	for (const auto& ca: cell_array){
		ofs << ca.size() << " ";
		for (int icell: ca){
			ofs << icell << " ";
		}
		ofs << std::endl;
	}

	ofs << "CELL_TYPES " << nc << std::endl;
	for (int v: cell_types){
		ofs << v << std::endl;
	}


	// data
	ofs << "POINT_DATA " << n_bases() << std::endl;
	for (auto& it: scalars){
		ofs << "SCALARS " << it.first << " double 1" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for (auto v: *it.second) ofs << v << std::endl;
	}
}
