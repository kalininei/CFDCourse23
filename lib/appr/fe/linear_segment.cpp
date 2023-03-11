#include "linear_segment.hpp"

LinearSegmentElement::LinearSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo), _len((coo[1] - coo[0]).x){}

std::vector<double> LinearSegmentElement::mass() const{
	return {_len/3, _len/6, _len/6, _len/3};
}

std::vector<double> LinearSegmentElement::stiff() const{
	double v = 1/_len;
	return {v, -v, -v, v};
}

std::vector<double> LinearSegmentElement::load() const{
	return {_len/2, _len/2};
}

LinearSegmentBoundaryElement::LinearSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo) {}

std::vector<double> LinearSegmentBoundaryElement::mass() const{
	return {1.0};
}
