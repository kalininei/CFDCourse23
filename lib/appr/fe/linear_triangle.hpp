#ifndef LINEAR_TRIANGLE_HPP
#define LINEAR_TRIANGLE_HPP

#include "appr/fe/aelement.hpp"
#include "appr/fe/linear_segment.hpp"

class LinearTriangleElement: public AElement{
public:
	LinearTriangleElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;
	std::vector<double> load() const override;
private:
	const double _modj;
};

class LinearTriangleBoundaryElement: public LinearSegmentElement{
public:
	LinearTriangleBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo, int direction);
};


#endif

