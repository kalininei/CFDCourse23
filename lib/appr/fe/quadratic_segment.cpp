#include "quadratic_segment.hpp"

namespace{
const SegmentQuad2 quad2;
const SegmentQuad3 quad3;
const SegmentQuad4 quad4;
}

QuadraticSegmentElement::QuadraticSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: ANumElement(nbases, coo, &quad2){}


std::array<double, 9> QuadraticSegmentElement::jacobi_matrix(Point p) const{
	return {
		point(1).x/2, 0, 0,
		0, 1, 0,
		0, 0, 1};
}

std::vector<double> QuadraticSegmentElement::bases(Point p) const{
	double xi_2 = p.x/2;
	double xi2_2 = p.x*p.x/2;
	return {
		-xi_2 + xi2_2,
		xi_2 + xi2_2,
		1 - 2*xi2_2};
}

std::vector<Point> QuadraticSegmentElement::grad_bases(Point p) const{
	double xi = p.x;
	return {
		-0.5 + xi,
		0.5 + xi,
		-2*xi};
}

QuadraticSegmentBoundaryElement::QuadraticSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo) {}

std::vector<double> QuadraticSegmentBoundaryElement::mass() const{
	return {1};
}
