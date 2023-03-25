#ifndef CUBIC_SEGMENT_HPP
#define CUBIC_SEGMENT_HPP

#include "appr/fe/anum_element.hpp"

class CubicSegmentElement: public ANumElement{
public:
	CubicSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
private:
	std::array<double, 9> jacobi_matrix(Point p) const override;
	std::vector<double> bases(Point p) const override;
	std::vector<Point> grad_bases(Point p) const override;
};


class CubicSegmentBoundaryElement: public AElement{
public:
	CubicSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
	std::vector<double> mass() const override;
};

#endif
