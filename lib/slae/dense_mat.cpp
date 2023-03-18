#include "dense_mat.hpp"
#include <cmath>

std::vector<double> dense_from_lower(const std::vector<double>& lower){
	int N = int((-1 + sqrt(8*lower.size() + 1))/2);
	std::vector<double> ret(N*N, 0);

	int k = 0;
	for (int j=0; j<N; ++j){
		for (int i=0; i<j; ++i){
			ret[N*j + i] = lower[k];
			ret[N*i + j] = lower[k];
			k++;
		}
		ret[N*j + j] = lower[k];
		k++;
	}
	return ret;
}

double determinant3(const std::array<double, 9>& mat){
	return mat[0]*(mat[4]*mat[8]-mat[5]*mat[7])
		-mat[1]*(mat[3]*mat[8]-mat[5]*mat[6])
		+mat[2]*(mat[3]*mat[7]-mat[4]*mat[6]);
}

std::array<double, 9> det_transpose_invert3(const std::array<double, 9>& a){
	return {
		a[4]*a[8]-a[5]*a[7], a[5]*a[6]-a[3]*a[8], a[3]*a[7]-a[4]*a[6],
		a[2]*a[7]-a[1]*a[8], a[0]*a[8]-a[2]*a[6], a[1]*a[6]-a[0]*a[7],
		a[1]*a[5]-a[2]*a[4], a[2]*a[3]-a[0]*a[5], a[0]*a[4]-a[1]*a[3]
	};
}

std::vector<Point> dense_mult_point(const std::array<double, 9>& mat, const std::vector<Point>& points){
	std::vector<Point> ret(points.size());

	for (size_t i=0; i<points.size(); ++i){
		const Point& p = points[i];
		ret[i].x = mat[0]*p.x + mat[1]*p.y + mat[2]*p.z;
		ret[i].y = mat[3]*p.x + mat[4]*p.y + mat[5]*p.z;
		ret[i].z = mat[6]*p.x + mat[7]*p.y + mat[8]*p.z;
	}

	return ret;
}

std::vector<double> outer_lower(const std::vector<double>& v1){
	int N = v1.size();
	std::vector<double> ret;
	for (int j=0; j<N; ++j){
		for (int i=0; i<j+1; ++i){
			ret.push_back(v1[i]*v1[j]);
		}
	}
	return ret;
}

std::vector<double> vector_dot_outer_lower(const std::vector<Vector>& v1){
	int N = v1.size();
	std::vector<double> ret;
	for (int j=0; j<N; ++j){
		for (int i=0; i<j+1; ++i){
			ret.push_back(vector_dot(v1[i], v1[j]));
		}
	}
	return ret;
}
