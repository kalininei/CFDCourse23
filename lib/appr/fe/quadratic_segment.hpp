#ifndef QUADRATIC_SEGMENT_HPP
#define QUADRATIC_SEGMENT_HPP

#include "appr/fe/anum_element.hpp"

class QuadraticSegmentElement: public ANumElement{
public:
	QuadraticSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
private:
	std::array<double, 9> jacobi_matrix(Point p) const override;
	std::vector<double> bases(Point p) const override;
	std::vector<Point> grad_bases(Point p) const override;
};


class QuadraticSegmentBoundaryElement: public AElement{
public:
	QuadraticSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
	std::vector<double> mass() const override;
};

#endif
