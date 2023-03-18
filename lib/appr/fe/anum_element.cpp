#include "anum_element.hpp"
#include "slae/dense_mat.hpp"

ANumElement::ANumElement(const std::vector<int>& nbases, const std::vector<Point>& coo, const IQuadratureRule* quad)
	: AElement(nbases, coo), _quad(quad){}

std::vector<double> ANumElement::mass() const{
	_THROW_NOT_IMP_;
}

std::vector<double> ANumElement::stiff() const{
	_THROW_NOT_IMP_;
}

std::vector<double> ANumElement::load() const{
	_THROW_NOT_IMP_;
}
