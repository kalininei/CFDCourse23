#ifndef BILINEAR_QUADRANGLE_HPP
#define BILINEAR_QUADRANGLE_HPP

#include "appr/fe/anum_element.hpp"
#include <vector>

class BilinearQuadrangleElement: public ANumElement{
public:
	BilinearQuadrangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
private:
	std::array<double, 9> jacobi_matrix(Point p) const override;
	std::vector<double> bases(Point p) const override;
	std::vector<Point> grad_bases(Point p) const override;
};

#endif

