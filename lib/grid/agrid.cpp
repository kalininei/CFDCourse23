#include "agrid.hpp"
#include <cassert>

AGrid::AGrid(int ndim): dim(ndim) { }


const AGridBoundary& AGrid::boundary(int ibnd) const{
	auto fnd = _boundary.find(ibnd);
	if (fnd == _boundary.end()) throw std::runtime_error(
		"Boundary " + std::to_string(ibnd) + " not found");
	return *(fnd->second);
}

std::vector<int> AGrid::btypes() const {
	std::vector<int> ret;
	for (auto bit: _boundary) ret.push_back(bit.first);
	return ret;
}

void AGrid::define_boundary(int btype, std::shared_ptr<AGridBoundary> boundary){
	_boundary[btype] = boundary;
}

void AGrid::define_boundary(int btype, std::function<bool(Point)> face_center_filter){
	// find boundary faces which are not assigned to any boundary
	std::vector<int> bfaces = boundary_faces();
	std::set<int> bfaces_set(bfaces.begin(), bfaces.end());
	for (auto it: _boundary){
		for (int iface: it.second->grid_face_indices()){
			bfaces_set.erase(iface);
		}
	}

	std::vector<int> bfaces_filtered(bfaces_set.begin(), bfaces_set.end());
	auto bnd = std::make_shared<AGridBoundary>(bfaces_filtered, face_center_filter, this);
	define_boundary(btype, bnd);
}

Point AGrid::face_center(int iface) const{
	int k = 0;
	Point ret {};
	for (int ipoint: tab_face_point(iface)){
		ret += point(ipoint);
		++k;
	}
	return ret / k;
}

Point AGrid::cell_center(int iface) const{
	int k = 0;
	Point ret {};
	for (int ipoint: tab_cell_point(iface)){
		ret += point(ipoint);
		++k;
	}
	return ret / k;
}

std::vector<int> AGrid::tab_cell_point(int icell) const{
	if (_cache.tab_cell_point.empty()){
		_cache.tab_cell_point.resize(n_cells());
		for (int icell=0; icell<n_cells(); ++icell){
			std::vector<int>& cp = _cache.tab_cell_point[icell];
			for (int iface: tab_cell_face(icell)){
				std::vector<int> fp = tab_face_point(iface);
				cp.insert(cp.end(), fp.begin(), fp.end());
			}

			std::sort(cp.begin(), cp.end());
			cp.resize(std::unique(cp.begin(), cp.end()) - cp.begin());
		}
	}
	return _cache.tab_cell_point[icell];
}

std::vector<int> AGrid::tab_cell_face(int icell) const{
	if (_cache.tab_cell_face.empty()){
		_cache.tab_cell_face.resize(n_cells());
		for (int iface=0; iface<n_faces(); ++iface){
			std::array<int, 2> fc = tab_face_cell(iface);
			if (fc[0] >= 0){
				_cache.tab_cell_face[fc[0]].push_back(iface);
			}
			if (fc[1] >= 0){
				_cache.tab_cell_face[fc[1]].push_back(iface);
			}
		}
	}
	return _cache.tab_cell_face[icell];
}

std::vector<int> AGrid::tab_cell_cell(int icell) const{
	if (_cache.tab_cell_cell.empty()){
		_cache.tab_cell_cell.resize(n_cells());
		for (int iface=0; iface<n_faces(); ++iface){
			std::array<int, 2> fc = tab_face_cell(iface);
			if (fc[0] >=0 && fc[1] >= 0){
				_cache.tab_cell_cell[fc[0]].push_back(fc[1]);
				_cache.tab_cell_cell[fc[1]].push_back(fc[0]);
			}
		}
	}

	return _cache.tab_cell_cell[icell];
}

std::vector<int> AGrid::tab_point_point(int ipoint) const{
	if (_cache.tab_point_point.empty()){
		std::vector<std::set<int>> pp(n_points());

		for (int icell=0; icell<n_cells(); ++icell){
			std::vector<int> cp = tab_cell_point(icell);
			for (int ipoint: cp){
				pp[ipoint].insert(cp.begin(), cp.end());
			}
		}

		for (int ipoint=0; ipoint<n_points(); ++ipoint){
			std::set<int>& s = pp[ipoint];
			s.erase(ipoint);
			_cache.tab_point_point.emplace_back(s.begin(), s.end());
		}
	}

	return _cache.tab_point_point[ipoint];
}

int AGrid::n_boundary_faces() const{
	if (_cache.n_boundary_faces < 0){
		_cache.n_boundary_faces = 0;
		for (int iface=0; iface < n_faces(); ++iface){
			std::array<int, 2> cc = tab_face_cell(iface);
			if (cc[0] < 0 || cc[1] < 0){
				_cache.n_boundary_faces += 1;
			}
		}
	}

	return _cache.n_boundary_faces;
}

int AGrid::n_internal_faces() const{
	if (_cache.n_internal_faces < 0){
		_cache.n_internal_faces = n_faces() - n_boundary_faces();
	}

	return _cache.n_internal_faces;
}

std::vector<int> AGrid::boundary_faces() const{
	if (_cache.boundary_faces.empty()){
		for (int iface=0; iface<n_faces(); ++iface){
			std::array<int, 2> cc = tab_face_cell(iface);
			if (cc[0] < 0 || cc[1] < 0){
				_cache.boundary_faces.push_back(iface);
			}
		}
	}
	return _cache.boundary_faces;
}

std::vector<int> AGrid::internal_faces() const{
	if (_cache.internal_faces.empty()){
		for (int iface=0; iface<n_faces(); ++iface){
			std::array<int, 2> cc = tab_face_cell(iface);
			if (cc[0] >= 0 && cc[1] >= 0){
				_cache.internal_faces.push_back(iface);
			}
		}
	}
	return _cache.internal_faces;
}

std::vector<int> AGrid::face_btypes() const{
	// initialize with -1 as internal faces
	std::vector<int> ret(n_faces(), -1);

	// zero to all boundary faces
	for (int iface: boundary_faces()){
		ret[iface] = 0;
	}

	// btype to all defined btypes
	for (auto bit: _boundary){
		for (int iface: bit.second->grid_face_indices()){
			ret[iface] = bit.first;
		}
	}

	return ret;
}

namespace{
void rearrange_index_sequence(std::vector<std::vector<int>>& ind){
	if (ind.size() < 2){
		return;
	}

	for (size_t i=0; i<ind.size()-1; ++i){
		int p0 = ind[i].back();
		for (size_t j=i+1; j<ind.size(); ++j){
			if (p0 == ind[j].front()){
				std::swap(ind[i+1], ind[j]);
				break;
			} else if (p0 == ind[j].back()){
				std::reverse(ind[j].begin(), ind[j].end());
				std::swap(ind[i+1], ind[j]);
				break;
			}
			if (j == ind.size() - 1){
				std::runtime_error("Failed to rearrange index sequence");
			}
		}
	}
}
}

std::vector<std::vector<int>> AGrid::vtk_cell_array() const{
	std::vector<std::vector<int>> ret;
	switch (dim){
		case 1:
		{
			for (int icell=0; icell<n_cells(); ++icell){
				ret.push_back(tab_cell_point(icell));
			}
			break;
		}
		case 2:
		{
			for (int icell=0; icell<n_cells(); ++icell){
				std::vector<std::vector<int>> cell_face_point;
				std::vector<int> faces = tab_cell_face(icell);
				for (int iface: faces){
					cell_face_point.push_back(tab_face_point(iface));
				}
				rearrange_index_sequence(cell_face_point);
				std::vector<int> p0;
				for (const auto& v: cell_face_point) p0.push_back(v[0]);
				if (tab_face_cell(faces[0])[0] != icell){
					std::reverse(p0.begin() + 1, p0.end());
				}
				ret.push_back(p0);
			}
			break;
		}
		case 3:
		{
			for (int icell=0; icell<n_cells(); ++icell){
				ret.emplace_back();
				auto& cret = ret.back();
				std::vector<int> faces = tab_cell_face(icell);
				cret.push_back(faces.size());
				for (int iface: faces){
					std::vector<int> points = tab_face_point(iface);
					cret.push_back(int(points.size()));
					cret.insert(cret.end(), points.begin(), points.end());
				}
			}
		}
	}

	return ret;
}

std::vector<int> AGrid::vtk_cell_types() const{
	switch (dim){
		case 1: return std::vector<int>(n_cells(), 3);
		case 2: return std::vector<int>(n_cells(), 7);
		case 3: return std::vector<int>(n_cells(), 42);
	}
	_THROW_UNREACHABLE_;
}

std::vector<std::vector<int>> AGrid::vtk_boundary_face_array() const{
	std::vector<std::vector<int>> ret;
	for (int iface: boundary_faces()){
		ret.push_back(tab_face_point(iface));
	}
	return ret;
}

std::vector<int> AGrid::vtk_boundary_face_types() const{
	switch (dim){
		case 1: return std::vector<int>(n_boundary_faces(), 1);
		case 2: return std::vector<int>(n_boundary_faces(), 3);
		case 3: return std::vector<int>(n_boundary_faces(), 7);
	}
	_THROW_UNREACHABLE_;
}

int AGrid::find_cell_index(Point p) const{
	_THROW_NOT_IMP_;
}

int AGrid::get_boundary_index_for_face(int iface) const{
	if (_cache.grid_to_bnd_face.size() == 0){
		std::vector<int> bfaces = boundary_faces();
		for (int i=0; i<(int)bfaces.size(); ++i){
			_cache.grid_to_bnd_face[bfaces[i]] = i;
		}
	}
	auto fnd = _cache.grid_to_bnd_face.find(iface);
	if (fnd != _cache.grid_to_bnd_face.end()){
		return fnd->second;
	} else {
		throw std::runtime_error("Face " + std::to_string(iface) + " is not boundary");
	}
}

std::array<double, 4> AGrid::face_plane(int iface) const{
	double A, B, C, D;
	Vector normal = face_normal(iface);
	Point pk = face_center(iface);
	A = normal.x;
	B = normal.y;
	C = normal.z;
	D = -(A*pk.x + B*pk.y + C*pk.z);
	return {A, B, C, D};
}
