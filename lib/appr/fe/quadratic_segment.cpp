#include "quadratic_segment.hpp"

namespace{
const SegmentQuad2 quad2;
const SegmentQuad3 quad3;
const SegmentQuad4 quad4;
}

QuadraticSegmentElement::QuadraticSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: ANumElement(nbases, coo, &quad2){}


std::array<double, 9> QuadraticSegmentElement::jacobi_matrix(Point p) const{
	_THROW_NOT_IMP_;
}

std::vector<double> QuadraticSegmentElement::bases(Point p) const{
	_THROW_NOT_IMP_;
}

std::vector<Point> QuadraticSegmentElement::grad_bases(Point p) const{
	_THROW_NOT_IMP_;
}

QuadraticSegmentBoundaryElement::QuadraticSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo) {}

std::vector<double> QuadraticSegmentBoundaryElement::mass() const{
	return {1};
}
