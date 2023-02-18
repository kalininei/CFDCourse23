#ifndef GEOM_HPP
#define GEOM_HPP

#include "common.hpp"

struct Point{
	double x=0, y=0, z=0;

	Point& operator+=(const Point& p){
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}

	Point& operator-=(const Point& p){
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}

	Point& operator*=(double c){
		x *= c;
		y *= c;
		z *= c;
		return *this;
	}

	Point& operator/=(double c){
		x /= c;
		y /= c;
		z /= c;
		return *this;
	}
};
using Vector = Point;

inline Point operator+(const Point& p1, const Point& p2){
	return {p1.x+p2.x, p1.y+p2.y, p1.z+p2.z};
}

inline Point operator-(const Point& p1, const Point& p2){
	return {p1.x-p2.x, p1.y-p2.y, p1.z-p2.z};
}

inline Point operator*(double c, const Point& p){
	return {c*p.x, c*p.y, c*p.z};
}

inline Point operator/(const Point& p, double c){
	return {p.x/c, p.y/c, p.z/c};
}

inline std::ostream& operator<<(std::ostream& s, const Point& p){
	s << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
	return s;
}

inline double vector_dot(const Vector& v1, const Vector& v2){
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline Vector vector_cross(const Vector& v1, const Vector& v2){
	return {
		v1.y*v2.z - v1.z*v2.y,
		-v1.x*v2.z + v1.z*v2.x,
		v1.x*v2.y - v1.y*v2.x
	};
}

inline double vector_dot_triple(const Vector& v1, const Vector& v2, const Vector& v3){
	return vector_dot(v1, vector_cross(v2, v3));
}

inline double vector_len(const Vector& v){
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline double vector_cos(const Vector& v1, const Vector& v2){
	return vector_dot(v1, v2)/vector_len(v1)/vector_len(v2);
}

inline double point_plane_distance(const Point& p, const Point& plane_point, const Vector& plane_normal){
	return vector_dot(plane_normal, p - plane_point);
}

#endif
