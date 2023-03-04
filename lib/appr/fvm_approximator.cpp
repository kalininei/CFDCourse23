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

double FvmApproximator::cell_to_face_center_normal_distance(int icell, int iface) const{
	// calculate plane equation
	std::array<double, 4> plane = _grid->face_plane(iface);
	double denum = std::sqrt(plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);
	Point pi = _grid->cell_center(icell);
	return std::abs(plane[0]*pi.x + plane[1]*pi.y + plane[2]*pi.z + plane[3])/denum;
}

std::vector<double> FvmApproximator::stiff() const{
	const CsrStencil& s = stencil();
	std::vector<double> ret(s.n_nonzero(), 0);

	// internal faces loop
	for (int iface: _grid->internal_faces()){
		int i, j;
		double area;
		double hik, hjk;

		std::array<int, 2> cells = _grid->tab_face_cell(iface);
		i = cells[0];
		j = cells[1];
		area = _grid->face_area(iface);
		
		// hik
		hik = cell_to_face_center_normal_distance(i, iface);
		hjk = cell_to_face_center_normal_distance(j, iface);
		
		double val = 1.0/(hik + hjk) * area;
		s.add_value(i, i, val, ret);
		s.add_value(j, j, val, ret);
		s.add_value(i, j, -val, ret);
		s.add_value(j, i, -val, ret);
	}

	// boundary faces loop
	int unk_k = _grid->n_cells();
	for (int iface: _grid->boundary_faces()){
		int i;
		double hik;
		double area;

		// i - basis for the adjacent cell
		std::array<int, 2> cells = _grid->tab_face_cell(iface);
		i = std::max(cells[0], cells[1]);

		// area
		area = _grid->face_area(iface);
		
		// hik
		hik = cell_to_face_center_normal_distance(i, iface);

		double val = 1.0/hik * area;

		s.add_value(i, i, val, ret);
		s.add_value(unk_k, unk_k, val, ret);
		s.add_value(i, unk_k, -val, ret);
		s.add_value(unk_k, i, -val, ret);

		++unk_k;
	}

	return ret;
}

std::map<int, std::vector<std::pair<int, Point>>> FvmApproximator::_build_boundary_bases() const{
	std::map<int, std::vector<std::pair<int, Point>>> ret;
	for (int btype: _grid->btypes()){
		std::vector<std::pair<int, Point>>& vec = ret[btype];
		const AGridBoundary& bnd = _grid->boundary(btype);
		for (int iface: bnd.grid_face_indices()){
			// index of the face amoung boundary faces
			int bnd_face_index = _grid->get_boundary_index_for_face(iface);
			// index of the basis
			int basis_index = _grid->n_cells() + bnd_face_index;
			vec.push_back({basis_index, _grid->face_center(iface)});
		}
	}
	return ret;
}

std::vector<double> FvmApproximator::_build_load_vector() const{
	std::vector<double> ret(n_bases(), 0);
	for (int i=0; i<_grid->n_cells(); ++i){
		ret[i] = _grid->cell_volume(i);
	}
	return ret;
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

double FvmApproximator::calculate_dudn(int btype, const std::vector<double>& v) const{
	double ret = 0;

	const AGridBoundary& bnd = _grid->boundary(btype);
	std::vector<int> face_indices = bnd.grid_face_indices();

	for (int iface: face_indices){
		// i - basis for face
		int bnd_index = _grid->get_boundary_index_for_face(iface);
		int i = _grid->n_cells() + bnd_index;

		// j - basis for adjacent cell
		std::array<int, 2> cells = _grid->tab_face_cell(iface);
		int icell = std::max(cells[0], cells[1]);

		// hij
		double hij = cell_to_face_center_normal_distance(icell, iface);

		double uj = v[icell];
		double ui = v[i];
		double area = _grid->face_area(iface);
		double val = (uj - ui)/hij*area;
		ret -= val;
	}

	return ret;
}

void FvmApproximator::apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const{
	const AGridBoundary& bnd = _grid->boundary(ibnd);
	std::vector<int> ifaces = bnd.grid_face_indices();

	for (int iface: ifaces){
		int bnd_index = _grid->get_boundary_index_for_face(iface);
		int ibasis = _grid->n_cells() + bnd_index;

		double q = q_func(_grid->face_center(iface));
		double area = _grid->face_area(iface);

		rhs[ibasis] -= q*area;
	};
}

void FvmApproximator::apply_bc_robin_to_stiff_lhs(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::vector<double>& stiff) const{

	const AGridBoundary& bnd = _grid->boundary(ibnd);
	std::vector<int> ifaces = bnd.grid_face_indices();
	const CsrStencil& s = stencil();

	for (int iface: ifaces){
		int bnd_index = _grid->get_boundary_index_for_face(iface);
		int ibasis = _grid->n_cells() + bnd_index;

		double alpha = alpha_func(_grid->face_center(iface));
		double area = _grid->face_area(iface);
		
		int diag_index = s.addr_index(ibasis, ibasis);
		stiff[diag_index] += area * alpha;
	};
}

void FvmApproximator::apply_bc_robin_to_stiff_rhs(
		int ibnd,
		std::function<double(Point)> beta_func,
		std::vector<double>& rhs) const{

	const AGridBoundary& bnd = _grid->boundary(ibnd);
	std::vector<int> ifaces = bnd.grid_face_indices();

	for (int iface: ifaces){
		int bnd_index = _grid->get_boundary_index_for_face(iface);
		int ibasis = _grid->n_cells() + bnd_index;

		double beta = beta_func(_grid->face_center(iface));
		double area = _grid->face_area(iface);
		
		rhs[ibasis] += area * beta;
	};
}

void FvmApproximator::apply_point_source(Point p, double flowrate, std::vector<double>& rhs) const{
	int ibasis = _grid->find_cell_index(p);
	if (ibasis < 0){
		throw std::runtime_error(
			"Failed to find cell for point: "
			+ std::to_string(p.x) + ", "
			+ std::to_string(p.y) + ", "
			+ std::to_string(p.z));
	}
	rhs[ibasis] += flowrate;
}
