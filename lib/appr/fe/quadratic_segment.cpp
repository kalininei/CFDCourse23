#include "quadratic_segment.hpp"

namespace{
const SegmentQuad2 quad2;
const SegmentQuad3 quad3;
const SegmentQuad4 quad4;
}

QuadraticSegmentElement::QuadraticSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: ANumElement(nbases, coo, &quad2){}


std::array<double, 9> QuadraticSegmentElement::jacobi_matrix(Point p) const{
	double len = point(1).x;
	return {
		len/2, 0, 0,
		0,     1, 0,
		0,     0, 1
	};
}

std::vector<double> QuadraticSegmentElement::bases(Point p) const{
	double xi = p.x;
	return {
		(-xi + xi*xi)/2,
		(xi + xi*xi)/2,
		1 - xi*xi
	};
}

std::vector<Point> QuadraticSegmentElement::grad_bases(Point p) const{
	double xi = p.x;
	return{
		Point(-0.5 + xi, 0, 0),
		Point(0.5 + xi, 0, 0),
		Point(-2*xi, 0, 0)
	};
}

QuadraticSegmentBoundaryElement::QuadraticSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo) {}

std::vector<double> QuadraticSegmentBoundaryElement::mass() const{
	return {1};
}
