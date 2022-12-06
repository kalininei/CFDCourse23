#include <numeric>
#include <fstream>
#include "appr/fvm_approximator.hpp"

int FvmApproximator::n_bases() const{
	return _grid->n_cells() + _grid->n_boundary_faces();
}

std::vector<double> FvmApproximator::approximate(std::function<double(Point)> fun) const{
	std::vector<double> ret(n_bases(), 0);
	int k = 0;
	for (int icell=0; icell < _grid->n_cells(); ++icell){
		ret[k++] = fun(_grid->cell_center(icell));
	}
	for (int iface: _grid->boundary_faces()){
		ret[k++] = fun(_grid->face_center(iface));
	}
	return ret;
}

std::shared_ptr<FvmApproximator> FvmApproximator::build(std::shared_ptr<AGrid> grid){
	return std::shared_ptr<FvmApproximator>(new FvmApproximator(grid));
}

FvmApproximator::FvmApproximator(std::shared_ptr<AGrid> grid): ASpatialApproximator(), _grid(grid){}

CsrStencil FvmApproximator::_build_stencil() const{
	std::vector<std::set<int>> sset(n_bases());

	// diagonals
	for (size_t i=0; i<sset.size(); ++i){
		sset[i].insert(i);
	}

	// cell->cell
	for (int icell=0; icell < _grid->n_cells(); ++icell)
	for (int icell2: _grid->tab_cell_cell(icell))
	if (icell2 > icell){
		sset[icell].insert(icell2);
		sset[icell2].insert(icell);
	}

	// boundary->cell
	int kface = 0;
	for (int iface: _grid->boundary_faces()){
		int face_basis = kface + _grid->n_cells();
		std::array<int, 2> icells = _grid->tab_face_cell(iface);
		int cell_basis = icells[0] >= 0 ? icells[0] : icells[1];

		sset[cell_basis].insert(face_basis);
		sset[face_basis].insert(cell_basis);

		++kface;
	}

	return CsrStencil::build(sset);
}

std::vector<double> FvmApproximator::mass() const{
	std::vector<double> ret(stencil().n_nonzero());
	std::vector<int> diags = stencil().addr_diag();

	for (int i=0; i<_grid->n_cells(); ++i){
		ret[diags[i]] = _grid->cell_volume(i);
	}
	for (int i=_grid->n_cells(); i<(int)diags.size(); ++i){
		ret[diags[i]] = 0.0;
	}

	return ret;
}

std::vector<double> FvmApproximator::stiff() const{
	std::vector<double> ret(stencil().n_nonzero(), 0);

	for (int iface: _grid->internal_faces()){
		std::array<int, 2> icells = _grid->tab_face_cell(iface);
		double face_area = _grid->face_area(iface);
		double h1 = point_plane_distance(_grid->cell_center(icells[0]), _grid->face_center(iface), _grid->face_normal(iface));
		double h2 = point_plane_distance(_grid->cell_center(icells[1]), _grid->face_center(iface), _grid->face_normal(iface));
		double coef = face_area / (std::abs(h1) + std::abs(h2));

		stencil().add_value(icells[0], icells[0], coef, ret);
		stencil().add_value(icells[1], icells[1], coef, ret);

		stencil().add_value(icells[0], icells[1], -coef, ret);
		stencil().add_value(icells[1], icells[0], -coef, ret);
	}

	int nbasis = _grid->n_cells();
	for (int iface: _grid->boundary_faces()){
		std::array<int, 2> icells = _grid->tab_face_cell(iface);
		int icell = (icells[0] >= 0) ? icells[0] : icells[1];
		double face_area = _grid->face_area(iface);
		double h = point_plane_distance(_grid->cell_center(icell), _grid->face_center(iface), _grid->face_normal(iface));
		double coef = face_area / std::abs(h);

		stencil().add_value(icell, icell, coef, ret);
		stencil().add_value(nbasis, nbasis, coef, ret);

		stencil().add_value(icell, nbasis, -coef, ret);
		stencil().add_value(nbasis, icell, -coef, ret);

		nbasis += 1;
	}

	return ret;
}

std::map<int, std::vector<std::pair<int, Point>>> FvmApproximator::_build_boundary_bases() const{
	std::map<int, std::vector<std::pair<int, Point>>> ret;

	std::vector<int> face_to_basis(_grid->n_faces(), -1);
	int it = _grid->n_cells();
	for (int iface: _grid->boundary_faces()){
		face_to_basis[iface] = it++;
	}
	
	for (int btype: _grid->btypes()){
		std::vector<std::pair<int, Point>> basis_info;
		
		for (int iface: _grid->boundary(btype).grid_face_indices()){
			int ibasis = face_to_basis[iface];
			Point fcenter = _grid->face_center(iface);
			basis_info.emplace_back(ibasis, fcenter);
		}

		ret[btype] = basis_info;
	}

	return ret;
}

std::vector<double> FvmApproximator::_build_load_vector() const{
	return mass();
}

void FvmApproximator::vtk_save_scalars(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const{
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
	ofs << "CELL_DATA " << nc << std::endl;
	for (auto& it: scalars){
		ofs << "SCALARS " << it.first << " double 1" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for (auto v: *it.second) ofs << v << std::endl;
	}
}
