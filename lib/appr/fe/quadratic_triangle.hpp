#ifndef QUADRATIC_TRIANGLE_HPP
#define QUADRATIC_TRIANGLE_HPP

#include "appr/fe/anum_element.hpp"
#include "appr/fe/quadratic_segment.hpp"

class QuadraticTriangleElement: public ANumElement{
public:
	QuadraticTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
private:
	std::array<double, 9> jacobi_matrix(Point p) const override;
	std::vector<double> bases(Point p) const override;
	std::vector<Point> grad_bases(Point p) const override;
};


class QuadraticTriangleBoundaryElement: public QuadraticSegmentElement{
public:
	QuadraticTriangleBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo, int direction);
	std::vector<double> mass() const override;
};

#endif
