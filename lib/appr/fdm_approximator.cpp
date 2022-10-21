#include <fstream>
#include "appr/fdm_approximator.hpp"

std::shared_ptr<FdmApproximator> FdmApproximator::build(std::shared_ptr<RegularGrid> grid){
	return std::shared_ptr<FdmApproximator>(new FdmApproximator(grid));
}

FdmApproximator::FdmApproximator(std::shared_ptr<RegularGrid> grid): ASpatialApproximator(), _grid(grid){}

int FdmApproximator::n_bases() const{
	return _grid->n_points();
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

void FdmApproximator::apply_bc_dirichlet_to_stiff_mat(int ibnd, std::vector<double>& stiff) const{
	const CsrStencil& s = stencil();
	const RegularGridBoundary& boundary = _grid->boundary(ibnd);

	for (int bp: boundary.point_indices()){
		s.set_unit_diagonal(bp, stiff);
	}
}

void FdmApproximator::apply_bc_dirichlet_to_stiff_vec(int ibnd, std::function<double(Point)> func, std::vector<double>& vec) const{
	const RegularGridBoundary& boundary = _grid->boundary(ibnd);

	for (int bp: boundary.point_indices()){
		vec[bp] = func(_grid->point(bp));
	}
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
	ofs << "z_COORDINATES " << _grid->nz() << " double" << std::endl;
	for (double c: _grid->zcoo()) ofs << c << std::endl;;

	// data
	ofs << "POINT_DATA " << _grid->n_points() << std::endl;
	for (auto& it: scalars){
		ofs << "SCALARS " << it.first << " double 1" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for (auto v: *it.second) ofs << v << std::endl;
	}


}
