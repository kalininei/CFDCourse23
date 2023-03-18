#include <fstream>
#include "grid/regular_grid.hpp"

namespace{

std::vector<double> uniform_range(double x0, double x1, int n){
	if (n == 0){
		throw std::runtime_error("Number of points in each direction should be >= 1");
	}
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
	if (ny > 1 && nx == 1){
		throw std::runtime_error("1D grids should by defined in x axis");
	}
	if (nz > 1 && (nx == 1 || ny == 1)){
		throw std::runtime_error("3D grids should be defined in xy axes");
	}

	if (ny == 1){
		return 1;
	} else if (nz == 1){
		return 2;
	} else {
		return 3;
	}
};

}

std::shared_ptr<ARegularGrid> ARegularGrid::build(
		int nx, double len_x,
		int ny, double len_y,
		int nz, double len_z){
	if (ny == 1){
		return std::make_shared<RegularGrid1>(nx, len_x);
	} else if (nz == 1){
		return std::make_shared<RegularGrid2>(nx, len_x, ny, len_y);
	} else {
		return std::make_shared<RegularGrid3>(nx, len_x, ny, len_y, nz, len_z);
	}
}

ARegularGrid::ARegularGrid(const std::vector<double>& xcoo,
                           const std::vector<double>& ycoo,
                           const std::vector<double>& zcoo)
		: AGrid(find_dim(xcoo.size(), ycoo.size(), zcoo.size())),
		  _xcoo(xcoo),
		  _ycoo(ycoo),
		  _zcoo(zcoo){ }

void ARegularGrid::define_reg_boundary(int btype, DirectionCode dircode){
	define_reg_boundary(btype, dircode, [](Point){ return true; });
}

void ARegularGrid::define_reg_boundary(int btype, DirectionCode dircode, std::function<bool(Point)> filter){
	define_boundary(btype, std::make_shared<RegularGridBoundary>(dircode, filter, this));
}

const RegularGridBoundary& ARegularGrid::reg_boundary(int ibnd) const{
	try{
		return dynamic_cast<const RegularGridBoundary&>(AGrid::boundary(ibnd));
	} catch (std::bad_cast){
		throw std::runtime_error(std::to_string(ibnd) + " boundary is not of RegularGridBoundary class");
	}
}

int ARegularGrid::n_points() const{
	return nx()*ny()*nz();
}

int ARegularGrid::nx() const{
	return _xcoo.size();
}

int ARegularGrid::ny() const{
	return _ycoo.size();
}

int ARegularGrid::nz() const{
	return _zcoo.size();
}

Point ARegularGrid::point(int point_index) const{
	std::array<int, 3> ijk = point_glob_to_ijk(point_index);
	return point_ijk(ijk[0], ijk[1], ijk[2]);
}

Point ARegularGrid::point_ijk(int ix, int iy, int iz) const{
	if (ix < 0 || ix >= nx()) throw std::runtime_error("X index is out of range: " + std::to_string(ix));
	if (iy < 0 || iy >= ny()) throw std::runtime_error("Y index is out of range: " + std::to_string(iy));
	if (iz < 0 || iz >= nz()) throw std::runtime_error("Z index is out of range: " + std::to_string(iz));
	return {_xcoo[ix], _ycoo[iy], _zcoo[iz]};
}

int ARegularGrid::point_ijk_to_glob(int ix, int iy, int iz) const{
	return ix + iy*nx() + iz*nx()*ny();
}

std::array<int, 3> ARegularGrid::point_glob_to_ijk(int point_index) const{
	std::array<int, 3> ret;

	int nxy = nx()*ny();

	ret[2] = point_index/nxy;
	ret[1] = (point_index - ret[2]*nxy) / nx();
	ret[0] = (point_index - ret[2]*nxy - ret[1]*nx());

	return ret;
}

int ARegularGrid::cell_ijk_to_glob(int ix, int iy, int iz) const{
	return ix + iy*(nx()-1) + iz*(nx()-1)*(ny()-1);
}

std::array<int, 3> ARegularGrid::cell_glob_to_ijk(int icell) const{
	std::array<int, 3> ret;

	ret[2] = icell / ((nx()-1) * (ny()-1));
	icell -= ret[2] * (nx()-1) * (ny()-1);
	ret[1] = icell / (nx()-1);
	ret[0] = icell % (nx()-1);

	return ret;
}

double ARegularGrid::hx(int ix) const{
	return _xcoo[ix+1] - _xcoo[ix];
}

double ARegularGrid::hy(int iy) const{
	return _ycoo[iy+1] - _ycoo[iy];
}

double ARegularGrid::hz(int iz) const{
	return _zcoo[iz+1] - _zcoo[iz];
}

RegularGrid1::RegularGrid1(int nx, double len_x)
	: RegularGrid1(uniform_range(0, len_x, nx)){
}

RegularGrid1::RegularGrid1(const std::vector<double>& xcoo): ARegularGrid(xcoo, {0}, {0}){
	if (dim != 1){
		throw std::runtime_error("Invalid RegularGrid1");
	}
}

RegularGrid2::RegularGrid2(int nx, double len_x, int ny, double len_y)
	: RegularGrid2(uniform_range(0, len_x, nx), uniform_range(0, len_y, ny)){
}

RegularGrid2::RegularGrid2(const std::vector<double>& xcoo, const std::vector<double>& ycoo)
		: ARegularGrid(xcoo, ycoo, {0}){
	if (dim != 2){
		throw std::runtime_error("Invalid RegularGrid2");
	}
}

RegularGrid3::RegularGrid3(
		int nx, double len_x,
		int ny, double len_y,
		int nz, double len_z)
		: RegularGrid3(uniform_range(0, len_x, nx), uniform_range(0, len_y, ny), uniform_range(0, len_z, nz)){
}

RegularGrid3::RegularGrid3(
		const std::vector<double>& xcoo,
		const std::vector<double>& ycoo,
		const std::vector<double>& zcoo)
		: ARegularGrid(xcoo, ycoo, zcoo){
	if (dim != 3){
		throw std::runtime_error("Invalid RegularGrid3");
	}
}

int RegularGrid1::face_ijk_to_glob(int ix, int iy, int iz, int idir) const{
	return ix;
}

std::array<int, 4> RegularGrid1::face_glob_to_ijk(int iface) const{
	return {iface, 0, 0, 0};
}

int RegularGrid2::face_ijk_to_glob(int ix, int iy, int iz, int idir) const{
	const int nx_xy = nx()*(ny()-1);
	switch (idir){
		case 0: return ix + iy*nx();
		case 1: return ix + iy*(nx()-1) + nx_xy;
	};
	_THROW_UNREACHABLE_;
}

std::array<int, 4> RegularGrid2::face_glob_to_ijk(int iface) const{
	const int nx_xy = nx()*(ny()-1);
	if (iface < nx_xy){
		return {iface % nx(), iface/nx(), 0, 0};
	} else {
		iface -= nx_xy;
		return {iface % (nx()-1), iface / (nx()-1), 0, 1};
	}
}


int RegularGrid3::face_ijk_to_glob(int ix, int iy, int iz, int idir) const{
	const int nx_xyz = nx()*(ny()-1)*(nz()-1);
	const int ny_xyz = (nx()-1)*ny()*(nz()-1);

	switch (idir){
		case 0: return ix + iy*nx() + iz*nx()*(ny()-1);
		case 1: return ix + iy*(nx()-1) + iz*(nx()-1)*ny() + nx_xyz;
		case 2: return ix + iy*(nx()-1) + iz*(nx()-1)*(ny()-1) + nx_xyz + ny_xyz;
	}
	_THROW_UNREACHABLE_;
}

std::array<int, 4> RegularGrid3::face_glob_to_ijk(int iface) const{
	const int nx_xyz = nx()*(ny()-1)*(nz()-1);
	const int ny_xyz = (nx()-1)*ny()*(nz()-1);

	if (iface < nx_xyz){
		int iz = iface / (nx()*(ny()-1));
		iface -= iz*nx()*(ny()-1);
		int iy = iface / nx();
		iface -= iy*nx();
		return {iface, iy, iz, 0};
	} else if (iface < nx_xyz + ny_xyz){
		iface -= nx_xyz;
		int iz = iface / ((nx()-1)*ny());
		iface -= iz * (nx()-1) * ny();
		int iy = iface / (nx()-1);
		iface -= iy * (nx()-1);
		return {iface, iy, iz, 1};
	} else {
		iface -= (nx_xyz + ny_xyz);
		int iz = iface / ((nx()-1)*(ny()-1));
		iface -= iz * (nx()-1)*(ny()-1);
		int iy = iface / (nx()-1);
		iface -= iy * (nx()-1);
		return {iface, iy, iz, 2};
	}
}

std::vector<int> RegularGrid1::tab_cell_point(int icell) const {
	return {icell, icell+1};
}

std::vector<int> RegularGrid2::tab_cell_point(int icell) const {
	std::array<int, 3> ijk = cell_glob_to_ijk(icell);
	return {
		point_ijk_to_glob(ijk[0], ijk[1], 0),
		point_ijk_to_glob(ijk[0]+1, ijk[1], 0),
		point_ijk_to_glob(ijk[0]+1, ijk[1]+1, 0),
		point_ijk_to_glob(ijk[0], ijk[1]+1, 0)};
}

std::vector<int> RegularGrid3::tab_cell_point(int icell) const {
	std::array<int, 3> ijk = cell_glob_to_ijk(icell);
	return {
		point_ijk_to_glob(ijk[0], ijk[1], ijk[2]),
		point_ijk_to_glob(ijk[0]+1, ijk[1], ijk[2]),
		point_ijk_to_glob(ijk[0]+1, ijk[1]+1, ijk[2]),
		point_ijk_to_glob(ijk[0], ijk[1]+1, ijk[2]),
		point_ijk_to_glob(ijk[0], ijk[1], ijk[2]+1),
		point_ijk_to_glob(ijk[0]+1, ijk[1], ijk[2]+1),
		point_ijk_to_glob(ijk[0]+1, ijk[1]+1, ijk[2]+1),
		point_ijk_to_glob(ijk[0], ijk[1]+1, ijk[2]+1)};
}


std::vector<int> RegularGrid1::tab_face_point(int iface) const {
	return {iface};
}

std::vector<int> RegularGrid2::tab_face_point(int iface) const {
	std::array<int, 4> locfac = face_glob_to_ijk(iface);
	int ipoint0 = point_ijk_to_glob(locfac[0], locfac[1], locfac[2]);
	switch (locfac[3]){
		case 0: return {ipoint0, ipoint0 + nx()};
		case 1: return {ipoint0+1, ipoint0};
	}
	_THROW_UNREACHABLE_;
}

std::vector<int> RegularGrid3::tab_face_point(int iface) const {
	std::array<int, 4> locfac = face_glob_to_ijk(iface);
	int ipoint0 = point_ijk_to_glob(locfac[0], locfac[1], locfac[2]);
	switch (locfac[3]){
		case 0: return {ipoint0, ipoint0 + nx(), ipoint0 + nx()*(1 + ny()), ipoint0 + nx()*ny()};
		case 1: return {ipoint0,  ipoint0 + nx()*ny(), ipoint0 + nx()*ny() + 1, ipoint0+1};
		case 2: return {ipoint0, ipoint0+1, ipoint0+1+nx(), ipoint0+nx()};
	}
	_THROW_UNREACHABLE_;
}

std::array<int, 2> RegularGrid1::tab_face_cell(int iface) const{
	if (iface == 0){
		return {-1, iface};
	} else if (iface == nx()-1){
		return {iface-1, -1};
	} else {
		return {iface-1, iface};
	}
}

std::array<int, 2> RegularGrid2::tab_face_cell(int iface) const{
	std::array<int, 4> faddr = face_glob_to_ijk(iface);
	switch (faddr[3]){
		case 0:
			return {faddr[0] == 0 ? -1 : cell_ijk_to_glob(faddr[0]-1, faddr[1], 0),
			        faddr[0] == nx()-1 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1], 0)};
		case 1:
			return {faddr[1] == 0 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1]-1, 0),
			        faddr[1] == ny()-1 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1], 0)};
	}
	_THROW_UNREACHABLE_;
}

std::array<int, 2> RegularGrid3::tab_face_cell(int iface) const{
	std::array<int, 4> faddr = face_glob_to_ijk(iface);
	switch (faddr[3]){
		case 0:
			return {faddr[0] == 0 ? -1 : cell_ijk_to_glob(faddr[0]-1, faddr[1], faddr[2]),
			        faddr[0] == nx()-1 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1], faddr[2])};
		case 1:
			return {faddr[1] == 0 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1]-1, faddr[2]),
			        faddr[1] == ny()-1 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1], faddr[2])};
		case 2:
			return {faddr[2] == 0 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1], faddr[2]-1),
			        faddr[2] == nz()-1 ? -1 : cell_ijk_to_glob(faddr[0], faddr[1], faddr[2])};
	}
	_THROW_UNREACHABLE_;
}

int RegularGrid1::n_boundary_faces() const{
	return 2;
}

int RegularGrid2::n_boundary_faces() const{
	return 2*(nx() + ny() - 2);
}

int RegularGrid3::n_boundary_faces() const{
	return 2*((nx()-1)*(ny()-1) + (nx()-1)*(nz()-1) + (ny()-1)*(nz()-1));
}

int RegularGrid1::n_cells() const{
	return nx()-1;
}

int RegularGrid2::n_cells() const{
	return (nx()-1)*(ny()-1);
}

int RegularGrid3::n_cells() const{
	return (nx()-1)*(ny()-1)*(nz()-1);
}

double RegularGrid1::cell_volume(int icell) const {
	return hx(icell);
}

double RegularGrid2::cell_volume(int icell) const {
	return hx(icell % nx()) * hy(icell / nx());
}

double RegularGrid3::cell_volume(int icell) const {
	std::array<int, 3> ijk = cell_glob_to_ijk(icell);
	return hx(ijk[0]) * hy(ijk[1]) * hz(ijk[2]);
}

double RegularGrid1::face_area(int iface) const {
	return 1;
}

double RegularGrid2::face_area(int iface) const {
	std::array<int, 4> faddr = face_glob_to_ijk(iface);
	if (faddr[3] == 0){
		return hy(faddr[1]);
	} else {
		return hx(faddr[0]);
	}
}

double RegularGrid3::face_area(int iface) const {
	std::array<int, 4> faddr = face_glob_to_ijk(iface);
	if (faddr[3] == 0){
		return hy(faddr[1])*hz(faddr[2]);
	} else if (faddr[3] == 1){
		return hx(faddr[0])*hz(faddr[2]);
	} else {
		return hy(faddr[1])*hz(faddr[2]);
	}
}

int RegularGrid1::n_faces() const{
	return nx();
}

int RegularGrid2::n_faces() const{
	return (nx()-1)*ny() + nx()*(ny()-1);
}

int RegularGrid3::n_faces() const{
	return nx() * (ny()-1) * (nz()-1) + (nx()-1) * ny() * (nz()-1) + (nx()-1) * (ny()-1) * nz();
}


Vector RegularGrid1::face_normal(int iface) const{
	return {1, 0, 0};
}

Vector RegularGrid2::face_normal(int iface) const{
	if (face_glob_to_ijk(iface)[3] == 0){
		return {1, 0, 0};
	} else {
		return {0, 1, 0};
	}
}

Vector RegularGrid3::face_normal(int iface) const{
	switch (face_glob_to_ijk(iface)[3]){
		case 0: return {1, 0, 0};
		case 1: return {0, 1, 0};
		case 2: return {0, 0, 1};
	}
	_THROW_UNREACHABLE_;
}
