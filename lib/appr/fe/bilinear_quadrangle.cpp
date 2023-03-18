#include "bilinear_quadrangle.hpp"
#include "fe_quadrature.hpp"
#include "slae/dense_mat.hpp"

namespace{
const SquareQuad2 quad2;
const SquareQuad3 quad3;
};

BilinearQuadrangleElement::BilinearQuadrangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
		: ANumElement(nbases, coo, &quad2){
	// check ordering
	double ar = vector_cross(coo[1] - coo[0], coo[2] - coo[0]).z;
	if (ar <= 0){
		throw std::runtime_error("Invalid ordering");
	}
}

std::array<double, 9> BilinearQuadrangleElement::jacobi_matrix(Point p) const{
	const double& xi = p.x;
	const double& eta = p.y;
	double j11 = (1-eta)/4 * point(1).x + (1+eta)/4 * point(2).x - (1+eta)/4 * point(3).x;
	double j12 = -(1+xi)/4 * point(1).x + (1+xi)/4 * point(2).x + (1-xi)/4 * point(3).x;
	double j21 = (1-eta)/4 * point(1).y + (1+eta)/4 * point(2).y - (1+eta)/4 * point(3).y;
	double j22 = -(1+xi)/4 * point(1).y + (1+xi)/4 * point(2).y + (1-xi)/4 * point(3).y;
	return {j11, j12, 0, j21, j22, 0, 0, 0, 1};
}

std::vector<double> BilinearQuadrangleElement::bases(Point p) const{
	const double& xi = p.x;
	const double& eta = p.y;
	return {(1-xi)*(1-eta)/4,
	        (1+xi)*(1-eta)/4,
	        (1+xi)*(1+eta)/4,
	        (1-xi)*(1+eta)/4};
}

std::vector<Point> BilinearQuadrangleElement::grad_bases(Point p) const{
	const double& xi = p.x;
	const double& eta = p.y;
	return {{-(1-eta)/4, -(1-xi)/4, 0 },
	        { (1-eta)/4, -(1+xi)/4, 0 },
	        { (1+eta)/4,  (1+xi)/4, 0 },
	        {-(1+eta)/4,  (1-xi)/4, 0 }};
}
