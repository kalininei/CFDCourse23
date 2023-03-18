#ifndef LINEAR_BUBLE_TRIANGLE_HPP
#define LINEAR_BUBLE_TRIANGLE_HPP

#include "appr/fe/anum_element.hpp"

class LinearBubleTriangleElement: public ANumElement{
public:
	LinearBubleTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
private:
	const double _modj;
	std::array<double, 9> jacobi_matrix(Point p) const override;
	std::vector<double> bases(Point p) const override;
	std::vector<Point> grad_bases(Point p) const override;
};

#endif
