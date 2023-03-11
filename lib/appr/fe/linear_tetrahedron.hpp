#ifndef LINEAR_TETRAHEDRON_HPP
#define LINEAR_TETRAHEDRON_HPP

#include "appr/fe/aelement.hpp"
#include "appr/fe/linear_triangle.hpp"

class LinearTetrahedronElement: public AElement{
public:
	LinearTetrahedronElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;
	std::vector<double> load() const override;
private:
	const double _modj;
};

class LinearTetrahedronBoundaryElement: public LinearTriangleElement{
public:
	LinearTetrahedronBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo, int direction);
};

#endif

