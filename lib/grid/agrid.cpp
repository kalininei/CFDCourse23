#include "agrid.hpp"

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

