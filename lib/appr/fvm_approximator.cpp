#include <numeric>
#include <fstream>
#include "appr/fvm_approximator.hpp"

int FvmApproximator::n_bases() const{
	return _grid->n_cells() + _grid->n_boundary_faces();
}

std::vector<double> FvmApproximator::approximate(std::function<double(Point)> fun) const{
	std::vector<double> ret(n_bases());

	int k=0;
	for (int icell=0; icell<_grid->n_cells(); ++icell){
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
	std::vector<std::set<int>> vec_set(n_bases());

	for (size_t i=0; i<vec_set.size(); ++i){
		vec_set[i].insert(i);
	}
	int k = 0;
	for (int icell=0; icell<_grid->n_cells(); ++icell){
		std::vector<int> cell_cell = _grid->tab_cell_cell(icell);
		std::set<int>& row_set = vec_set[k++];
		for (int cc: cell_cell) row_set.insert(cc);
	}
	for (int iface: _grid->boundary_faces()){
		std::array<int, 2> face_cell = _grid->tab_face_cell(iface);
		int icell = face_cell[0];
		if (icell < 0) icell = face_cell[1];
		int irow_face = k++;
		vec_set[icell].insert(irow_face);
		vec_set[irow_face].insert(icell);
		std::cout << icell << " " << irow_face << std::endl;
	}
	return CsrStencil::build(vec_set);
}

std::vector<double> FvmApproximator::mass() const{
	const CsrStencil& s = stencil();
	std::vector<double> ret(s.n_nonzero(), 0);

	for (int icell=0; icell<_grid->n_cells(); ++icell){
		int addr = s.addr_index(icell, icell);
		ret[addr] = _grid->cell_volume(icell);
	}

	return ret;
}

std::vector<double> FvmApproximator::stiff() const{
	_THROW_NOT_IMP_;
}

std::map<int, std::vector<std::pair<int, Point>>> FvmApproximator::_build_boundary_bases() const{
	_THROW_NOT_IMP_;
}

std::vector<double> FvmApproximator::_build_load_vector() const{
	_THROW_NOT_IMP_;
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
