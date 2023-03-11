#ifndef LINEAR_SEGMENT_HPP
#define LINEAR_SEGMENT_HPP

#include "appr/fe/aelement.hpp"

class LinearSegmentElement: public AElement{
public:
	LinearSegmentElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;
	std::vector<double> load() const override;
private:
	const double _len;
};


class LinearSegmentBoundaryElement: public AElement{
public:
	LinearSegmentBoundaryElement(const std::vector<int>& nbases, const std::vector<Point>& coo);
	std::vector<double> mass() const override;
};

#endif
