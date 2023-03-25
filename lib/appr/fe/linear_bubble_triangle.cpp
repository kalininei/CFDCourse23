#include "linear_bubble_triangle.hpp"

namespace{
const TriQuad2 quad2;
const TriQuad3 quad3;
const TriQuad4 quad4;
}

LinearBubbleTriangleElement::LinearBubbleTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo):
	ANumElement(nbases, coo, &quad4){}
	
std::array<double, 9> LinearBubbleTriangleElement::jacobi_matrix(Point p) const{
	Point p1 = point(1);
	Point p2 = point(2);
	return {
		p1.x, p2.x, 0,
		p1.y, p2.y, 0,
		0,    0,    1};
}

std::vector<double> LinearBubbleTriangleElement::bases(Point p) const{
	double xi = p.x;
	double eta = p.y;
	double v = 9*xi*eta*(1-xi-eta);
	return {
		1 - xi - eta - v,
		xi - v,
		eta - v,
		3*v};
}

std::vector<Point> LinearBubbleTriangleElement::grad_bases(Point p) const{
	double xi = p.x;
	double eta = p.y;
	double dvd1 = 9*eta - 18*xi*eta - 9*eta*eta;
	double dvd2 = 9*xi - 18*xi*eta - 9*xi*xi;
	return {
		Point(-1 - dvd1, -1-dvd2, 0),
		Point(1-dvd1, -dvd2, 0),
		Point(-dvd1, 1-dvd2, 0),
		Point(3*dvd1, 3*dvd2, 0)
	};
}
