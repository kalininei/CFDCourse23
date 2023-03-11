#include "aelement.hpp"

namespace{
std::vector<Point> shift_point_vector(const std::vector<Point>& coo){
	std::vector<Point> ret;
	ret.push_back({0, 0, 0});

	for (size_t i=1; i<coo.size(); ++i){
		ret.push_back(coo[i] - coo[0]);
	}

	return ret;
}
}

AElement::AElement(const std::vector<int>& bases, const std::vector<Point>& coo)
	: _bases(bases), _coo(shift_point_vector(coo)), _point0(coo[0])
{
}

// total number of local bases
int AElement::n_bases() const{
	return (int)_bases.size();
}

// global index of local basis
int AElement::i_basis(int ilocal) const{
	return _bases[ilocal];
}

// total number of points in the cell
int AElement::n_points() const{
	return (int)_coo.size();
}

// shifted point coordinate. First is always 0
Point AElement::point(int ilocal) const{
	return _coo[ilocal];
}

std::vector<double> AElement::mass() const{
	_THROW_NOT_IMP_;
}

std::vector<double> AElement::stiff() const{
	_THROW_NOT_IMP_;
}

std::vector<double> AElement::load() const{
	_THROW_NOT_IMP_;
}
