#include "appr/spatial_approximator.hpp"
#include <fstream>
#include <numeric>

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

const std::vector<double>& ASpatialApproximator::load_vector() const{
	if (_cache.load_vector.empty()){
		_cache.load_vector = _build_load_vector();
	}
	return _cache.load_vector;
}

void ASpatialApproximator::apply_bc_dirichlet(
		int ibnd,
		std::function<double(Point)> val_func,
		std::vector<double>& lhs, std::vector<double>& rhs) const{
	const CsrStencil& s = stencil();

	const std::vector<std::pair<int, Point>>& bnd_points = boundary_bases(ibnd);

	if (!lhs.empty())
	for (std::pair<int, Point> it: bnd_points){
		s.set_unit_diagonal(it.first, lhs);
	}

	if (!rhs.empty())
	for (std::pair<int, Point> it: bnd_points){
		rhs[it.first] = val_func(it.second);
	}
}

void ASpatialApproximator::apply_bc_dirichlet_lhs(int ibnd, std::vector<double>& lhs) const{
	std::vector<double> rhs;
	return apply_bc_dirichlet(ibnd, [](Point){ return 0; }, lhs, rhs);
}

void ASpatialApproximator::apply_bc_dirichlet_rhs(int ibnd, std::function<double(Point)> val_func, std::vector<double>& rhs) const{
	std::vector<double> lhs;
	return apply_bc_dirichlet(ibnd, val_func, lhs, rhs);
}

const std::vector<std::pair<int, Point>>& ASpatialApproximator::boundary_bases(int ibnd) const{
	if (_cache.boundary_bases.empty()){
		_cache.boundary_bases = _build_boundary_bases();
	}
	return _cache.boundary_bases[ibnd];
}

void ASpatialApproximator::apply_bc_neumann_to_stiff(
		int ibnd,
		std::function<double(Point)> q_func,
		std::vector<double>& rhs) const{
	_THROW_NOT_IMP_;
}

void ASpatialApproximator::apply_bc_robin_to_stiff(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::function<double(Point)> beta_func,
		std::vector<double>& stiff, std::vector<double>& rhs) const{
	_THROW_NOT_IMP_;
}

double ASpatialApproximator::integrate(const std::vector<double>& v) const{
	const std::vector<double>& lv = load_vector();
	double ret = 0;

	for (int i=0; i<n_bases(); ++i){
		ret += v[i] * lv[i];
	}

	return ret;
}

double ASpatialApproximator::domain_volume() const{
	if (_cache.volume == 0){
		const std::vector<double>& lv = load_vector();
		_cache.volume = std::accumulate(lv.begin(), lv.end(), 0.0);
	}
	return _cache.volume;
}

double ASpatialApproximator::norm_max(const std::vector<double>& v) const{
	// only for nodal v
	if ((int)v.size() < n_bases()) throw std::invalid_argument("norm_max()");
	double ret = v[0];
	for (int i=1; i<n_bases(); ++i){
		ret = std::max(ret, std::abs(v[i]));
	}
	return ret;
}

double ASpatialApproximator::norm_2(const std::vector<double>& v) const{
	// !!! only for nodal v
	if ((int)v.size() < n_bases()) throw std::invalid_argument("norm_2()");

	std::vector<double> v2(n_bases());
	for (int i=0; i<n_bases(); ++i){
		v2[i] = v[i]*v[i];
	}
	return std::sqrt(integrate(v2));
}

double ASpatialApproximator::rms(const std::vector<double>& v) const{
	return 1.0/std::sqrt(domain_volume()) * norm_2(v);
}

double ASpatialApproximator::rms(const std::vector<double>& a, const std::vector<double>& b) const{
	std::vector<double> diff(a);
	for (size_t i=0; i<a.size(); ++i) diff[i] -= b[i];
	return rms(diff);
}

double ASpatialApproximator::rms(const std::vector<double>& a, std::function<double(Point)> b) const{
	return rms(a, approximate(b));
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
