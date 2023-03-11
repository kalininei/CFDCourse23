#include "linear_segment.hpp"

LinearSegmentElement::LinearSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo), _len((coo[1] - coo[0]).x){}

std::vector<double> LinearSegmentElement::mass() const{
	_THROW_NOT_IMP_;
}

std::vector<double> LinearSegmentElement::stiff() const{
	_THROW_NOT_IMP_;
}

std::vector<double> LinearSegmentElement::load() const{
	_THROW_NOT_IMP_;
}

LinearSegmentBoundaryElement::LinearSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo) {}

std::vector<double> LinearSegmentBoundaryElement::mass() const{
	return {1.0};
}
