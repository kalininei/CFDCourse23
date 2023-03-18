#include "anum_element.hpp"
#include "slae/dense_mat.hpp"

ANumElement::ANumElement(const std::vector<int>& nbases, const std::vector<Point>& coo, const IQuadratureRule* quad)
	: AElement(nbases, coo), _quad(quad){}

std::vector<double> ANumElement::mass() const{
	auto fun = [&](Point p)->std::vector<double>{
		std::array<double, 9> j = jacobi_matrix(p);
		std::vector<double> phi = bases(p);
		std::vector<double> ij = outer_lower(phi);
		double m = determinant3(j);
		for (double& v: ij) v*=m;
		return ij;
	};
	std::vector<double> mm = _quad->integrate(fun);
	return dense_from_lower(mm);
}

std::vector<double> ANumElement::stiff() const{
	auto fun = [&](Point p)->std::vector<double>{
		std::array<double, 9> j = jacobi_matrix(p);
		std::array<double, 9> jt = det_transpose_invert3(j);
		std::vector<Point> dphi_dxi = grad_bases(p);
		std::vector<Point> dphi_dx  = dense_mult_point(jt, dphi_dxi);
		std::vector<double> ij = vector_dot_outer_lower(dphi_dx);
		double m = determinant3(j);
		for (double& v: ij) v/=m;
		return ij;
	};
	std::vector<double> mm = _quad->integrate(fun);
	return dense_from_lower(mm);
}

std::vector<double> ANumElement::load() const{
	auto fun = [&](Point p)->std::vector<double>{
		std::array<double, 9> j = jacobi_matrix(p);
		double m = determinant3(j);
		std::vector<double> phi = bases(p);
		for (double& v: phi) v *= m;
		return phi;
	};
	return _quad->integrate(fun);
}
