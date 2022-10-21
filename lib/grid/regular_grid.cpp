#include <fstream>
#include "grid/regular_grid.hpp"

namespace{

std::vector<double> uniform_range(double x0, double x1, int n){
	if (n == 0) return {0};
	std::vector<double> ret;
	for (int i=0; i<n; ++i){
		ret.push_back(x0 + i*(x1-x0)/(n-1));
	}
	return ret;
};

int find_dim(int nx, int ny, int nz){
	if (nx < 1 || ny < 1 || nz < 1) {
		throw std::runtime_error("Regular grid coordinates dimensions should be >= 1");
	}

	int ret = 0;
	if (nx > 1) ret += 1;
	if (ny > 1) ret += 1;
	if (nz > 1) ret += 1;

	return ret;
};

}

RegularGrid::RegularGrid(
		double len_x, int nx,
		double len_y, int ny,
		double len_z, int nz):
	RegularGrid(uniform_range(0, len_x, nx),
	            uniform_range(0, len_y, ny),
	            uniform_range(0, len_z, nz)) { }



RegularGrid::RegularGrid(const std::vector<double>& xcoo,
                         const std::vector<double>& ycoo,
                         const std::vector<double>& zcoo):
		dim(find_dim(xcoo.size(), ycoo.size(), zcoo.size())),
		_xcoo(xcoo),
		_ycoo(ycoo),
		_zcoo(zcoo){ }

void RegularGrid::define_boundary(int btype, DirectionCode dircode){
	define_boundary(btype, dircode, [](Point){ return true; });
}

void RegularGrid::define_boundary(int btype, DirectionCode dircode, std::function<bool(Point)> filter){
	int istart, iend, jstart, jend, kstart, kend;
	switch (dircode){
		case DirectionCode::X_PLUS:
			istart = _xcoo.size()-1;
			iend = _xcoo.size();
			jstart = 0;
			jend = _ycoo.size();
			kstart = 0;
			kend = _zcoo.size();
			break;
		case DirectionCode::X_MINUS:
			istart = 0;
			iend = 1;
			jstart = 0;
			jend = _ycoo.size();
			kstart = 0;
			kend = _zcoo.size();
			break;
		case DirectionCode::Y_PLUS:
			istart = 0;
			iend = _xcoo.size();
			jstart = _ycoo.size()-1;
			jend = _ycoo.size();
			kstart = 0;
			kend = _zcoo.size();
			break;
		case DirectionCode::Y_MINUS:
			istart = 0;
			iend = _xcoo.size();
			jstart = 0;
			jend = 1;
			kstart = 0;
			kend = _zcoo.size();
			break;
		case DirectionCode::Z_PLUS:
			istart = 0;
			iend = _xcoo.size();
			jstart = 0;
			jend = _ycoo.size();
			kstart = _zcoo.size()-1;
			kend = _zcoo.size();
			break;
		case DirectionCode::Z_MINUS:
			istart = 0;
			iend = _xcoo.size();
			jstart = 0;
			jend = _ycoo.size();
			kstart = 0;
			kend = 1;
			break;
	}

	std::vector<int> pindices;
	for (int k=kstart; k<kend; ++k)
	for (int j=jstart; j<jend; ++j)
	for (int i=istart; i<iend; ++i){
		int g = ijk_to_glob(i, j, k);
		if (filter(point_ijk(i, j, k))){
			pindices.push_back(g);
		}
	}

	for (auto& bit: _boundary){
		RegularGridBoundary& b = *bit.second;
		if (b.direction == dircode){
			b.remove_points(pindices);
		}
	}

	_boundary[btype] = std::make_shared<RegularGridBoundary>(dircode, pindices);
}

int RegularGrid::n_points() const{
	return nx()*ny()*nz();
}

int RegularGrid::nx() const{
	return _xcoo.size();
}

int RegularGrid::ny() const{
	return _ycoo.size();
}

int RegularGrid::nz() const{
	return _zcoo.size();
}

Point RegularGrid::point(int point_index) const{
	std::array<int, 3> ijk = glob_to_ijk(point_index);
	return point_ijk(ijk[0], ijk[1], ijk[2]);
}

Point RegularGrid::point_ijk(int ix, int iy, int iz) const{
	if (ix < 0 || ix >= nx()) throw std::runtime_error("X index is out of range: " + std::to_string(ix));
	if (iy < 0 || iy >= ny()) throw std::runtime_error("Y index is out of range: " + std::to_string(iy));
	if (iz < 0 || iz >= nz()) throw std::runtime_error("Z index is out of range: " + std::to_string(iz));
	return {_xcoo[ix], _ycoo[iy], _zcoo[iz]};
}

int RegularGrid::ijk_to_glob(int ix, int iy, int iz) const{
	return ix + iy*nx() + iz*nx()*ny();
}

std::array<int, 3> RegularGrid::glob_to_ijk(int point_index) const{
	std::array<int, 3> ret;

	int nxy = nx()*ny();

	ret[2] = point_index/nxy;
	ret[1] = (point_index - ret[2]*nxy) / nx();
	ret[0] = (point_index - ret[2]*nxy - ret[1]*nx());

	return ret;
}

double RegularGrid::hx(int ix) const{
	return _xcoo[ix+1] - _xcoo[ix];
}

double RegularGrid::hy(int iy) const{
	return _ycoo[iy+1] - _ycoo[iy];
}

double RegularGrid::hz(int iz) const{
	return _zcoo[iz+1] - _zcoo[iz];
}

const RegularGridBoundary& RegularGrid::boundary(int ibnd) const{
	auto fnd = _boundary.find(ibnd);
	if (fnd == _boundary.end()) throw std::runtime_error(
		"Boundary " + std::to_string(ibnd) + " not found");
	return *(fnd->second);
}
