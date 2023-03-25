#include "quadratic_triangle.hpp"

namespace{
const TriQuad2 quad2;
const TriQuad3 quad3;
const TriQuad4 quad4;
}

QuadraticTriangleElement::QuadraticTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: ANumElement(nbases, coo, &quad2){}

namespace{

std::vector<int> build_bnd_nbases(const std::vector<int>& nbases, int direction){
	if (direction == 1){
		return nbases;
	} else {
		std::vector<int> ret(nbases);
		std::reverse(ret.begin(), ret.end()-1);
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

QuadraticTriangleBoundaryElement::QuadraticTriangleBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo, int direction)
	: QuadraticSegmentElement(build_bnd_nbases(nbases, direction), build_bnd_coo(coo, direction))
{}
