#include "linear_triangle.hpp"

namespace{

double compute_modj(const std::vector<Point>& coo){
	double j11 = coo[1].x - coo[0].x;
	double j12 = coo[2].x - coo[0].x;
	double j21 = coo[1].y - coo[0].y;
	double j22 = coo[2].y - coo[0].y;
	return j22*j11 - j12*j21;
}

};

LinearTriangleElement::LinearTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo), _modj(compute_modj(coo)){}

std::vector<double> LinearTriangleElement::mass() const{
	double v = _modj/24.0;
	return {2*v, v, v,
		v, 2*v, v,
		v, v, 2*v};
}

std::vector<double> LinearTriangleElement::stiff() const{
	std::vector<double> ret(9);

	double j11 = point(1).x;
	double j12 = point(2).x;
	double j21 = point(1).y;
	double j22 = point(2).y;

	ret[0] = ((j21 - j22)*(j21 - j22) + (j11 - j12)*(j11 - j12))/2/_modj;
	ret[1] = ( j22*(j21 - j22) + j12*(j11 - j12) )/2/_modj;
	ret[2] = (-j11*(j11 - j12)  - j21*(j21 - j22) )/2/_modj;
	ret[4] = ( j22*j22   + j12*j12  )/2/_modj;
	ret[5] = (-j21*j22   - j12*j11  )/2/_modj;
	ret[8] = ( j11*j11   + j21*j21  )/2/_modj;
	ret[3] = ret[1];
	ret[6] = ret[2];
	ret[7] = ret[5];

	return ret;

}

std::vector<double> LinearTriangleElement::load() const{
	double v = _modj/6;
	return {v, v, v};
}

namespace{

std::vector<int> build_bnd_nbases(const std::vector<int>& nbases, int direction){
	if (direction == 1){
		return nbases;
	} else {
		std::vector<int> ret(nbases);
		std::reverse(ret.begin(), ret.end());
		return ret;
	}
}

std::vector<Point> build_bnd_coo(const std::vector<Point>& coo, int direction){
	std::vector<Point> ret(2);
	ret[0].x = 0;
	ret[0].y = 0;
	ret[0].z = 0;

	ret[1].x = vector_len(coo[1] - coo[0]);
	ret[1].y = 0;
	ret[1].z = 0;

	return ret;
}

}

LinearTriangleBoundaryElement::LinearTriangleBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo, int direction)
	: LinearSegmentElement(build_bnd_nbases(nbases, direction), build_bnd_coo(coo, direction))
{}
