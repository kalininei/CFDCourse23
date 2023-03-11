#include "linear_tetrahedron.hpp"

namespace{

double compute_modj(const std::vector<Point>& coo){
	double j11 = coo[1].x - coo[0].x;
	double j12 = coo[2].x - coo[0].x;
	double j13 = coo[3].x - coo[0].x;
	double j21 = coo[1].y - coo[0].y;
	double j22 = coo[2].y - coo[0].y;
	double j23 = coo[3].y - coo[0].y;
	double j31 = coo[1].z - coo[0].z;
	double j32 = coo[2].z - coo[0].z;
	double j33 = coo[3].z - coo[0].z;

	return j11*(j22*j33 - j23*j32) - j12*(j21*j33 - j23*j31) + j13*(j21*j32 - j22*j31);
}

};

LinearTetrahedronElement::LinearTetrahedronElement(const std::vector<int>& nbases, const std::vector<Point>& coo)
	: AElement(nbases, coo), _modj(compute_modj(coo)){}

std::vector<double> LinearTetrahedronElement::mass() const{
	double v = _modj/120;
	return {2*v, v, v, v,
	        v, 2*v, v, v,
	        v, v, 2*v, v,
	        v, v, v, 2*v};
}

std::vector<double> LinearTetrahedronElement::load() const{
	double v = _modj/12;
	return {v, v, v, v};
}

std::vector<double> LinearTetrahedronElement::stiff() const{
	double j11 = point(1).x;
	double j12 = point(2).x;
	double j13 = point(3).x;
	double j21 = point(1).y;
	double j22 = point(2).y;
	double j23 = point(3).y;
	double j31 = point(1).z;
	double j32 = point(2).z;
	double j33 = point(3).z;

	double m11 = j22*j33-j23*j32;
	double m12 = j23*j31-j21*j33;
	double m13 = j21*j32-j22*j31;
	double m21 = j13*j32-j12*j33;
	double m22 = j11*j33-j13*j31;
	double m23 = j12*j31-j11*j32;
	double m31 = j12*j23-j13*j22;
	double m32 = j13*j21-j11*j23;
	double m33 = j11*j22-j12*j21;

	// dphi1/dx
	double d1d1 = m11*(-1) + m12*(-1) + m13*(-1);
	// dphi1/dy
	double d1d2 = m21*(-1) + m22*(-1) + m23*(-1);
	// dphi1/dz
	double d1d3 = m31*(-1) + m32*(-1) + m33*(-1);

	// dphi2/dx
	double d2d1 = m11*(1) + m12*(0) + m13*(0);
	// dphi2/dy
	double d2d2 = m21*(1) + m22*(0) + m23*(0);
	// dphi2/dz
	double d2d3 = m31*(1) + m32*(0) + m33*(0);

	double d3d1 = m11*(0) + m12*(1) + m13*(0);
	double d3d2 = m21*(0) + m22*(1) + m23*(0);
	double d3d3 = m31*(0) + m32*(1) + m33*(0);

	double d4d1 = m11*(0) + m12*(0) + m13*(1);
	double d4d2 = m21*(0) + m22*(0) + m23*(1);
	double d4d3 = m31*(0) + m32*(0) + m33*(1);

	std::vector<double> ret(16);
	ret[0] = d1d1*d1d1 + d1d2*d1d2 + d1d3*d1d3;
	ret[1] = d1d1*d2d1 + d1d2*d2d2 + d1d3*d2d3;
	ret[2] = d1d1*d3d1 + d1d2*d3d2 + d1d3*d3d3;
	ret[3] = d1d1*d4d1 + d1d2*d4d2 + d1d3*d4d3;

	ret[5] = d2d1*d2d1 + d2d2*d2d2 + d2d3*d2d3;
	ret[6] = d2d1*d3d1 + d2d2*d3d2 + d2d3*d3d3;
	ret[7] = d2d1*d4d1 + d2d2*d4d2 + d2d3*d4d3;

	ret[10] = d3d1*d3d1 + d3d2*d3d2 + d3d3*d3d3;
	ret[11] = d3d1*d4d1 + d3d2*d4d2 + d3d3*d4d3;

	ret[15] = d4d1*d4d1 + d4d2*d4d2 + d4d3*d4d3;

	// symmetrical
	ret[4] = ret[1];
	ret[8] = ret[2];
	ret[9] = ret[6];
	ret[12] = ret[3];
	ret[13] = ret[7];
	ret[14] = ret[11];

	for (auto& v: ret){
		v/=(6*_modj);
	}

	return ret;
}

namespace{

std::vector<int> build_bnd_nbases(std::vector<int> nbases, int direction){
	if (direction == 1){
		return nbases;
	} else {
		std::reverse(nbases.begin() + 1, nbases.end());
		return nbases;
	}
}

std::vector<Point> build_bnd_coo(std::vector<Point> coo, int direction){
	if (direction == -1){
		std::reverse(coo.begin() + 1, coo.end());
	}
	Vector x1 = coo[1] - coo[0];
	Vector ytmp = coo[2] - coo[0];
	Vector z1 = vector_cross(x1, ytmp);
	Vector y1 = vector_cross(z1, x1);
	x1 /= vector_len(x1);
	y1 /= vector_len(y1);

	std::vector<Point> ret(coo.size());
	ret[1].x = vector_dot(coo[1] - coo[0], x1);
	ret[1].y = vector_dot(coo[1] - coo[0], y1);
	ret[2].x = vector_dot(coo[2] - coo[0], x1);
	ret[2].y = vector_dot(coo[2] - coo[0], y1);

	return ret;
}

}

LinearTetrahedronBoundaryElement::LinearTetrahedronBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo, int direction)
	: LinearTriangleElement(build_bnd_nbases(nbases, direction), build_bnd_coo(coo, direction))
{
}
