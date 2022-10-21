#include "appr/spatial_approximator.hpp"
#include <fstream>

const CsrStencil& ASpatialApproximator::stencil() const{
	if (_cache.stencil.n_rows() == 0){
		_cache.stencil = _build_stencil();
	}
	return _cache.stencil;
}

std::vector<double> ASpatialApproximator::stiff() const{
	_THROW_NOT_IMP_;
}

std::vector<double> ASpatialApproximator::mass() const{
	_THROW_NOT_IMP_;
}

std::vector<double> ASpatialApproximator::transport(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const{
	_THROW_NOT_IMP_;
}

std::vector<double> ASpatialApproximator::transport_upwind(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const{
	_THROW_NOT_IMP_;
}

void ASpatialApproximator::apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff_mat) const{
	_THROW_NOT_IMP_;
}

void ASpatialApproximator::apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> dir_values,
                                                           std::vector<double>& stiff_vec) const{
	_THROW_NOT_IMP_;
}

void ASpatialApproximator::vtk_save_scalar(std::string filepath, const std::vector<double>& scalar, std::string datacap) const{
	return _vtk_save_scalar(filepath, {{datacap, &scalar}});
}

void ASpatialApproximator::vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const{
	return _vtk_save_scalar(filepath, scalars);
}

void ASpatialApproximator::_vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const{
	_THROW_NOT_IMP_;
}
