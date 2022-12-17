#include <numeric>
#include <fstream>
#include "appr/fvm_approximator.hpp"

int FvmApproximator::n_bases() const{
	_THROW_NOT_IMP_;
}

std::vector<double> FvmApproximator::approximate(std::function<double(Point)> fun) const{
	_THROW_NOT_IMP_;
}

std::shared_ptr<FvmApproximator> FvmApproximator::build(std::shared_ptr<AGrid> grid){
	_THROW_NOT_IMP_;
}

FvmApproximator::FvmApproximator(std::shared_ptr<AGrid> grid): ASpatialApproximator(), _grid(grid){}

CsrStencil FvmApproximator::_build_stencil() const{
	_THROW_NOT_IMP_;
}

std::vector<double> FvmApproximator::mass() const{
	_THROW_NOT_IMP_;
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
