#include <fstream>
#include "appr/fdm_approximator.hpp"

std::shared_ptr<FdmApproximator> FdmApproximator::build(std::shared_ptr<RegularGrid> grid){
	return std::shared_ptr<FdmApproximator>(new FdmApproximator(grid));
}

FdmApproximator::FdmApproximator(std::shared_ptr<RegularGrid> grid): ASpatialApproximator(), _grid(grid){}

int FdmApproximator::n_bases() const{
	return _grid->n_points();
}

std::map<int, std::vector<std::pair<int, Point>>> FdmApproximator::_build_boundary_bases() const{
	std::map<int, std::vector<std::pair<int, Point>>> ret;

	// sort btypes in reverse order to give higher priority to large btype values
	std::vector<int> btypes = _grid->btypes();
	std::sort(btypes.begin(), btypes.end());
	std::reverse(btypes.begin(), btypes.end());

	// using std::set to prevent doubling of points in different bc faces
	std::set<int> used_points;

	for (int btype: btypes){
		std::vector<std::pair<int, Point>>& vec = ret[btype];
		const RegularGridBoundary& bnd = _grid->boundary(btype);
		for (int ipoint: bnd.point_indices()){
			std::pair<std::set<int>::const_iterator, bool> ires = used_points.insert(ipoint);
			// if point was not used before
			if (ires.second){
				vec.push_back({ipoint, _grid->point(ipoint)});
			}
		}
	}

	return ret;
}

std::vector<double> FdmApproximator::approximate(std::function<double(Point)> func) const{
	std::vector<double> ret(n_bases());

	int k = 0;
	for (int iz=0; iz<_grid->nz(); ++iz)
	for (int iy=0; iy<_grid->ny(); ++iy)
	for (int ix=0; ix<_grid->nx(); ++ix){
		ret[k++] = func(_grid->point_ijk(ix, iy, iz));
	}

	return ret;
}

std::vector<double> FdmApproximator::stiff() const{
	const CsrStencil& s = stencil();
	std::vector<double> ret(s.n_nonzero());

	int nx = _grid->nx();
	int ny = _grid->ny();
	int nz = _grid->nz();
	int nxy = nx*ny;

	double hplus, hminus, vplus, vminus;
	int ijk;

	// x
	if (_grid->nx() > 1)
	for (int iz=0; iz<nz; ++iz)
	for (int iy=0; iy<ny; ++iy){
		//left
		ijk = _grid->ijk_to_glob(0, iy, iz);
		hplus = _grid->hx(0);
		vplus = 2/hplus/hplus;
		ret[s.addr_index(ijk, ijk)] += vplus;
		ret[s.addr_index(ijk, ijk+1)] -= vplus;

		// internal
		for (int ix=1; ix<nx-1; ++ix){
			ijk = _grid->ijk_to_glob(ix, iy, iz);
			hplus = _grid->hx(ix);
			hminus = _grid->hx(ix-1);
			vplus = 2/(hplus+hminus)/hplus;
			vminus = 2/(hplus+hminus)/hminus;
			ret[s.addr_index(ijk, ijk)] += vplus + vminus;
			ret[s.addr_index(ijk, ijk-1)] -= vminus;
			ret[s.addr_index(ijk, ijk+1)] -= vplus;
		}

		// right
		ijk = _grid->ijk_to_glob(nx-1, iy, iz);
		hminus = _grid->hx(nx-2);
		vminus = 2/hminus/hminus;
		ret[s.addr_index(ijk, ijk)] += vminus;
		ret[s.addr_index(ijk, ijk-1)] -= vminus;
	}

	// y
	if (_grid->ny() > 1)
	for (int iz=0; iz<nz; ++iz)
	for (int ix=0; ix<nx; ++ix){
		//left
		ijk = _grid->ijk_to_glob(ix, 0, iz);
		hplus = _grid->hy(0);
		vplus = 2/hplus/hplus;
		ret[s.addr_index(ijk, ijk)] += vplus;
		ret[s.addr_index(ijk, ijk+nx)] -= vplus;

		// internal
		for (int iy=1; iy<ny-1; ++iy){
			ijk = _grid->ijk_to_glob(ix, iy, iz);
			hplus = _grid->hy(iy);
			hminus = _grid->hy(iy-1);
			vplus = 2/(hplus+hminus)/hplus;
			vminus = 2/(hplus+hminus)/hminus;
			ret[s.addr_index(ijk, ijk)] += vplus + vminus;
			ret[s.addr_index(ijk, ijk-nx)] -= vminus;
			ret[s.addr_index(ijk, ijk+nx)] -= vplus;
		}

		// right
		ijk = _grid->ijk_to_glob(ix, ny-1, iz);
		hminus = _grid->hy(ny-2);
		vminus = 2/hminus/hminus;
		ret[s.addr_index(ijk, ijk)] += vminus;
		ret[s.addr_index(ijk, ijk-nx)] -= vminus;
	}

	// z
	if (_grid->nz() > 1)
	for (int iy=0; iy<ny; ++iy)
	for (int ix=0; ix<nx; ++ix){
		//left
		ijk = _grid->ijk_to_glob(ix, iy, 0);
		hplus = _grid->hz(0);
		vplus = 2/hplus/hplus;
		ret[s.addr_index(ijk, ijk)] += vplus;
		ret[s.addr_index(ijk, ijk+nxy)] -= vplus;

		// internal
		for (int iz=1; iz<nz-1; ++iz){
			ijk = _grid->ijk_to_glob(ix, iy, iz);
			hplus = _grid->hz(iz);
			hminus = _grid->hz(iz-1);
			vplus = 2/(hplus+hminus)/hplus;
			vminus = 2/(hplus+hminus)/hminus;
			ret[s.addr_index(ijk, ijk)] += vplus + vminus;
			ret[s.addr_index(ijk, ijk-nxy)] -= vminus;
			ret[s.addr_index(ijk, ijk+nxy)] -= vplus;
		}

		// right
		ijk = _grid->ijk_to_glob(ix, iy, nz-1);
		hminus = _grid->hz(nz-2);
		vminus = 2/hminus/hminus;
		ret[s.addr_index(ijk, ijk)] += vminus;
		ret[s.addr_index(ijk, ijk-nxy)] -= vminus;
	}

	return ret;
}

std::vector<double> FdmApproximator::mass() const{
	const CsrStencil& s = stencil();
	std::vector<double> ret(s.n_nonzero(), 0);

	for (int i=0; i<s.n_rows(); ++i){
		ret[s.addr()[i]] = 1.0;
	}

	return ret;
}

CsrStencil FdmApproximator::_build_stencil() const{
	std::vector<std::set<int>> ij(n_bases());

	int nx = _grid->nx();
	int ny = _grid->ny();
	int nz = _grid->nz();
	int nxy = nx*ny;

	for (int iz=0; iz<nz; ++iz)
	for (int iy=0; iy<ny; ++iy)
	for (int ix=0; ix<nx; ++ix){
		int ijk = ix + iy*nx + iz*nxy;
		std::set<int>& row = ij[ijk];
		// diagonal
		row.insert(ijk);
		// +x
		if (ix < nx-1) row.insert(ijk+1);
		// -x
		if (ix > 0) row.insert(ijk-1);
		// +y
		if (iy < ny-1) row.insert(ijk+nx);
		// -y
		if (iy > 0) row.insert(ijk-nx);
		// +z
		if (iz < nz-1) row.insert(ijk+nxy);
		// -z
		if (iz > 0) row.insert(ijk-nxy);
	}

	return CsrStencil::build(ij);
}

std::vector<double> FdmApproximator::transport(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const{
	_THROW_NOT_IMP_;
}

std::vector<double> FdmApproximator::transport_upwind(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const{
	const CsrStencil& s = stencil();
	int nx = _grid->nx();
	int nxy = _grid->nx() * _grid->ny();

	std::vector<double> ret(s.n_nonzero(), 0);

	for (int iz=0; iz<_grid->nz(); ++iz)
	for (int iy=0; iy<_grid->ny(); ++iy)
	for (int ix=0; ix<_grid->nx(); ++ix){
		int ipoint = _grid->ijk_to_glob(ix, iy, iz);

		// x direction
		if (vx[ipoint] > 0 && ix > 0){
			double v = vx[ipoint]/_grid->hx(ix-1);
			ret[s.addr_index(ipoint, ipoint)] += v;
			ret[s.addr_index(ipoint, ipoint-1)] -= v;
		} else if (vx[ipoint] < 0 && ix < _grid->nx()-1){
			double v = vx[ipoint]/_grid->hx(ix);
			ret[s.addr_index(ipoint, ipoint)] -= v;
			ret[s.addr_index(ipoint, ipoint+1)] += v;
		}
	
		// y direction
		if (vy[ipoint] > 0 && iy > 0){
			double v = vy[ipoint]/_grid->hy(iy-1);
			ret[s.addr_index(ipoint, ipoint)] += v;
			ret[s.addr_index(ipoint, ipoint-nx)] -= v;
		} else if (vy[ipoint] < 0 && iy < _grid->ny()-1){
			double v = vy[ipoint]/_grid->hy(iy);
			ret[s.addr_index(ipoint, ipoint)] -= v;
			ret[s.addr_index(ipoint, ipoint+nx)] += v;
		}

		// z direction
		if (vz[ipoint] > 0 && iz > 0){
			double v = vz[ipoint]/_grid->hz(iz-1);
			ret[s.addr_index(ipoint, ipoint)] += v;
			ret[s.addr_index(ipoint, ipoint-nxy)] -= v;
		} else if (vz[ipoint] < 0 && iz < _grid->nz()-1){
			double v = vz[ipoint]/_grid->hz(iz);
			ret[s.addr_index(ipoint, ipoint)] -= v;
			ret[s.addr_index(ipoint, ipoint+nxy)] += v;
		}
	}

	return ret;
}

std::vector<double> FdmApproximator::_build_load_vector() const{
	std::vector<double> ret(n_bases(), 0);

	for (int iz=0; iz<_grid->nz(); ++iz)
	for (int iy=0; iy<_grid->ny(); ++iy)
	for (int ix=0; ix<_grid->nx(); ++ix){
		double hx = 0;
		double hy = 0;
		double hz = 0;

		if (ix > 0) hx += _grid->hx(ix-1)/2;
		if (ix < _grid->nx()-1) hx += _grid->hx(ix)/2;
		if (iy > 0) hy += _grid->hy(iy-1)/2;
		if (iy < _grid->ny()-1) hy += _grid->hy(iy)/2;
		if (iz > 0) hz += _grid->hz(iz-1)/2;
		if (iz < _grid->nz()-1) hz += _grid->hz(iz)/2;

		if (hx == 0) hx = 1;
		if (hy == 0) hy = 1;
		if (hz == 0) hz = 1;

		int ipoint = _grid->ijk_to_glob(ix, iy, iz);
		ret[ipoint] = hx*hy*hz;
	}

	return ret;
};

void FdmApproximator::_vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const{
	std::ofstream ofs(filepath);

	// header
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "cfdlib" << std::endl;
	ofs << "ASCII" << std::endl;

	// grid
	ofs << "DATASET RECTILINEAR_GRID" << std::endl;
	ofs << "DIMENSIONS " << _grid->nx() << " " << _grid->ny() << " " << _grid->nz() << std::endl;
	ofs << "X_COORDINATES " << _grid->nx() << " double" << std::endl;
	for (double c: _grid->xcoo()) ofs << c << std::endl;;
	ofs << "Y_COORDINATES " << _grid->ny() << " double" << std::endl;
	for (double c: _grid->ycoo()) ofs << c << std::endl;;
	ofs << "Z_COORDINATES " << _grid->nz() << " double" << std::endl;
	for (double c: _grid->zcoo()) ofs << c << std::endl;;

	// data
	ofs << "POINT_DATA " << _grid->n_points() << std::endl;
	for (auto& it: scalars){
		ofs << "SCALARS " << it.first << " double 1" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for (auto v: *it.second) ofs << v << std::endl;
	}


}

double FdmApproximator::boundary_h(int ibnd) const{
	DirectionCode code = _grid->boundary(ibnd).direction;

	switch (code){
		case DirectionCode::X_MINUS: return _grid->hx(0);
		case DirectionCode::X_PLUS:  return _grid->hx(_grid->nx()-2);
		case DirectionCode::Y_MINUS: return _grid->hy(0);
		case DirectionCode::Y_PLUS:  return _grid->hy(_grid->ny()-2);
		case DirectionCode::Z_MINUS: return _grid->hz(0);
		case DirectionCode::Z_PLUS:  return _grid->hz(_grid->nz()-2);
	}

	_THROW_NOT_IMP_;
}

void FdmApproximator::apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const{
	// 1. получить узлы (inode), принадлежащие границе ibnd
	const std::vector<std::pair<int, Point>>& nodes = boundary_bases(ibnd);

	// 2. добавить в rhs[inode] необходимое слагаемое
	// 2.1 вычисление множителя h
	double h = boundary_h(ibnd);

	// 2.2 В цикле по узлам границы вычисляем du/dn и добавляем в rhs
	for (const auto& b: nodes){
		double dudn = -q_func(b.second);
		double value = 2/h * dudn;
		rhs[b.first] += value;
	}
}

void FdmApproximator::apply_bc_robin_to_stiff_lhs(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::vector<double>& stiff) const{
	// 1. получить узлы (inode), принадлежащие границе ibnd
	const std::vector<std::pair<int, Point>>& nodes = boundary_bases(ibnd);
	
	// 2. вычисление множителя h
	double h = boundary_h(ibnd);

	// 3. Подставить в диагональ матрицы
	for (const std::pair<int, Point>& ninfo: nodes){
		// compute value 2*alpha/h
		double value = 2*alpha_func(ninfo.second)/h;
		// insert into diagonal
		int diagonal_addr = stencil().addr_index(ninfo.first, ninfo.first);
		stiff[diagonal_addr] += value;
	}
}

void FdmApproximator::apply_bc_robin_to_stiff_rhs(
		int ibnd,
		std::function<double(Point)> beta_func,
		std::vector<double>& rhs) const{
	// 1. получить узлы (inode), принадлежащие границе ibnd
	const std::vector<std::pair<int, Point>>& nodes = boundary_bases(ibnd);
	
	// 2. вычисление множителя h
	double h = boundary_h(ibnd);

	// 3. Подставить в правую часть
	for (const std::pair<int, Point>& ninfo: nodes){
		// compute value 2*beta/h
		double value = 2*beta_func(ninfo.second)/h;
		rhs[ninfo.first] += value;
	}
}
