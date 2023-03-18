#include "anum_element.hpp"
#include "slae/dense_mat.hpp"

ANumElement::ANumElement(const std::vector<int>& nbases, const std::vector<Point>& coo, const IQuadratureRule* quad)
	: AElement(nbases, coo), _quad(quad){}

std::vector<double> ANumElement::mass() const{
	auto fun = [&](Point p){
		std::vector<double> phi_i = bases(p);
		std::array<double, 9> jac = jacobi_matrix(p);
		double modj = determinant3(jac);
		std::vector<double> ij = outer_lower(phi_i);
		for (double& v: ij) v *= modj;
		return ij;
	};
	std::vector<double> lower_mass = _quad->integrate(fun);
	return dense_from_lower(lower_mass);
}

std::vector<double> ANumElement::stiff() const{
	auto fun = [&](Point p){
		std::vector<Point> dphi = grad_bases(p);
		std::array<double, 9> jac = jacobi_matrix(p);
		std::array<double, 9> jac_ti = det_transpose_invert3(jac);
		std::vector<Point> dphi_dx = dense_mult_point(jac_ti, dphi);
		std::vector<double> ij = vector_dot_outer_lower(dphi_dx);
		double modj = determinant3(jac);
		for (double& v: ij) v /= modj;
		return ij;
	};
	std::vector<double> lower_stiff = _quad->integrate(fun);
	return dense_from_lower(lower_stiff);
}

std::vector<double> ANumElement::load() const{
	auto fun = [&](Point p){
		std::vector<double> phi_i = bases(p);
		std::array<double, 9> jac = jacobi_matrix(p);
		double modj = determinant3(jac);

		for (double& v: phi_i) v *= modj;

		return phi_i;
	};
	return _quad->integrate(fun);
}
