#include "cubic_segment.hpp"

namespace{
const SegmentQuad2 quad2;
const SegmentQuad3 quad3;
const SegmentQuad4 quad4;
}

CubicSegmentElement::CubicSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: ANumElement(nbases, coo, &quad3){}


std::array<double, 9> CubicSegmentElement::jacobi_matrix(Point p) const{
	return {
		point(1).x/2, 0, 0,
		0, 1, 0,
		0, 0, 1};
}

std::vector<double> CubicSegmentElement::bases(Point p) const{
	double xi = p.x;
	double xi2 = p.x*p.x;
	double xi3 = p.x*p.x*p.x;
	double m = 1.0/16.0;
	return {
		m*(-1 + xi + 9*xi2 - 9*xi3),
		m*(-1 - xi + 9*xi2 + 9*xi3),
		m*( 9 - 27*xi - 9*xi2 + 27*xi3),
		m*( 9 + 27*xi - 9*xi2 - 27*xi3)};
}

std::vector<Point> CubicSegmentElement::grad_bases(Point p) const{
	double xi = p.x;
	double xi2 = p.x*p.x;
	double m = 1.0/16.0;
	return {
		m*(1 + 18*xi - 27*xi2),
		m*(-1 +18*xi + 27*xi2),
		m*(-27 - 18*xi + 81*xi2),
		m*( 27 - 18*xi - 81*xi2)};
}

CubicSegmentBoundaryElement::CubicSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo) {}

std::vector<double> CubicSegmentBoundaryElement::mass() const{
	return {1};
}
