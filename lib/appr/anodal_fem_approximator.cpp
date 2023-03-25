#include <sstream>
#include <fstream>
#include <numeric>
#include "anodal_fem_approximator.hpp"
#include "debug.hpp"

ANodalFemApproximator::ANodalFemApproximator(std::shared_ptr<AGrid> grid): _grid(grid){
}

void ANodalFemApproximator::initialize(){
	for (int icell=0; icell < _grid->n_cells(); ++icell){
		_elements.emplace_back(build_element(icell));
	}

	for (int iface: _grid->boundary_faces()){
		_boundary_elements[iface] = std::shared_ptr<AElement>(build_boundary_element(iface));
	}
}

int ANodalFemApproximator::n_elements() const {
	return (int)_elements.size();
}

std::vector<double> ANodalFemApproximator::approximate(std::function<double(Point)> func) const{
	std::vector<double> ret(n_bases());

	for (int i=0; i<n_bases(); ++i){
		ret[i] = func(node(i));
	}

	return ret;
}

CsrStencil ANodalFemApproximator::_build_stencil() const{
	std::vector<std::set<int>> s(n_bases());

	for (std::shared_ptr<AElement> el: _elements){
		for (int i=0; i<el->n_bases(); ++i){
			int gi = el->i_basis(i);
			s[gi].insert(gi);

			for (int j=0; j<i; ++j){
				int gj = el->i_basis(j);
				s[gi].insert(gj);
				s[gj].insert(gi);
			}
		}
	}

	return CsrStencil::build(s);
}

std::vector<double> ANodalFemApproximator::mass() const{
	const CsrStencil& s = stencil();
	std::vector<double> ret(s.n_nonzero(), 0);

	for (std::shared_ptr<AElement> el: _elements){
		std::vector<double> local_mass = el->mass();
		for (int j=0; j<el->n_bases(); ++j){
			for (int i=0; i<el->n_bases(); ++i){
				double v = local_mass[j*el->n_bases() + i];
				int irow = el->i_basis(j);
				int icol = el->i_basis(i);
				int ind = s.addr_index(irow, icol);
				ret[ind] += v;
			}
		}
	}
	//dbg::print(s, ret);

	return ret;
}

std::vector<double> ANodalFemApproximator::stiff() const{
	const CsrStencil& s = stencil();
	std::vector<double> ret(s.n_nonzero(), 0);

	for (std::shared_ptr<AElement> el: _elements){
		std::vector<double> local_stiff = el->stiff();
		for (int j=0; j<el->n_bases(); ++j){
			for (int i=0; i<el->n_bases(); ++i){
				double v = local_stiff[j*el->n_bases() + i];
				int irow = el->i_basis(j);
				int icol = el->i_basis(i);
				int ind = s.addr_index(irow, icol);
				ret[ind] += v;
			}
		}
	}
	//dbg::print(s, ret);
	return ret;
}

std::vector<double> ANodalFemApproximator::_build_load_vector() const{
	std::vector<double> ret(n_bases(), 0);

	for (std::shared_ptr<AElement> el: _elements){
		std::vector<double> local_vec = el->load();
		for (int j=0; j<el->n_bases(); ++j){
			double v = local_vec[j];
			int ind = el->i_basis(j);
			ret[ind] += v;
		}
	}

	return ret;
}

void ANodalFemApproximator::apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const{
	const AGridBoundary& bnd = _grid->boundary(ibnd);
	std::vector<int> ifaces = bnd.grid_face_indices();
	for (int iface: ifaces){
		const std::shared_ptr<AElement>& el = _boundary_elements.find(iface)->second;
		std::vector<double> local_mass_mat = el->mass();
		std::vector<double> local_mass_vec(el->n_bases(), 0);
		for (int irow=0; irow<el->n_bases(); ++irow){
			for (int icol=0; icol<el->n_bases(); ++icol){
				double v_icol = q_func(node(el->i_basis(icol)));
				local_mass_vec[irow] -= local_mass_mat[irow*el->n_bases() + icol] * v_icol;
			}
		}

		add_local_vector(1, local_mass_vec, el.get(), rhs);
	}
}

void ANodalFemApproximator::apply_bc_robin_to_stiff_lhs(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::vector<double>& stiff) const {

	const AGridBoundary& bnd = _grid->boundary(ibnd);
	std::vector<int> ifaces = bnd.grid_face_indices();
	for (int iface: ifaces){
		const std::shared_ptr<AElement>& el = _boundary_elements.find(iface)->second;
		std::vector<double> local_mass_mat = el->mass();
		double alpha_e = alpha_func(_grid->face_center(iface));
		add_local_matrix(alpha_e, local_mass_mat, el.get(), stiff);
	}
}

void ANodalFemApproximator::apply_bc_robin_to_stiff_rhs(
		int ibnd,
		std::function<double(Point)> beta_func,
		std::vector<double>& rhs) const {

	auto minus_beta = [&](Point p){return -beta_func(p); };
	apply_bc_neumann_to_stiff(ibnd, minus_beta, rhs);
}


void ANodalFemApproximator::add_local_matrix(double coef, const std::vector<double>& lmat, const AElement* elem, std::vector<double>& gmat) const{
	const CsrStencil& s = stencil();
	for (int j=0; j<elem->n_bases(); ++j){
		int gj = elem->i_basis(j);
		for (int i=0; i<elem->n_bases(); ++i){
			int gi = elem->i_basis(i);
			double v = lmat[j*elem->n_bases() + i];
			int plain_index = s.addr_index(gj, gi);
			gmat[plain_index] += coef * v;
		}
	}
}

void ANodalFemApproximator::add_local_vector(double coef, const std::vector<double>& lvec, const AElement* elem, std::vector<double>& gvec) const{
	for (int i=0; i<elem->n_bases(); ++i){
		int gi = elem->i_basis(i);
		gvec[gi] += coef*lvec[i];
	}
}

std::map<int, std::vector<std::pair<int, Point>>> ANodalFemApproximator::_build_boundary_bases() const{
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
		for (int iface: bnd.grid_face_indices()){
			auto el = _boundary_elements.find(iface)->second;
			for (int ilocal=0; ilocal < el->n_bases(); ++ilocal){
				int iglobal = el->i_basis(ilocal);
				std::pair<std::set<int>::const_iterator, bool> ires = used_points.insert(iglobal);
				// if point was not used before
				if (ires.second){
					vec.push_back({iglobal, node(iglobal)});
				}
			}
		}
	}

	return ret;
}

void ANodalFemApproximator::vtk_save_scalars(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const{
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
	ofs << "POINT_DATA " << _grid->n_points() << std::endl;
	for (auto& it: scalars){
		ofs << "SCALARS " << it.first << " double 1" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for (int i = 0; i<_grid->n_points(); ++i){
			ofs << (*it.second)[i] << std::endl;
		}
	}
}
