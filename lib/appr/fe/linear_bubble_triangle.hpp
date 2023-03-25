#ifndef LINEAR_BUBBLE_TRIANGLE_HPP
#define LINEAR_BUBBLE_TRIANGLE_HPP

#include "appr/fe/anum_element.hpp"

class LinearBubbleTriangleElement: public ANumElement{
public:
	LinearBubbleTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
private:
	std::array<double, 9> jacobi_matrix(Point p) const override;
	std::vector<double> bases(Point p) const override;
	std::vector<Point> grad_bases(Point p) const override;
};

#endif
