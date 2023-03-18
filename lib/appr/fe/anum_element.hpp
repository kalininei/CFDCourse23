#ifndef _ANUM_ELEMENT_HPP
#define _ANUM_ELEMENT_HPP

#include "appr/fe/aelement.hpp"
#include "appr/fe/fe_quadrature.hpp"

// element with numerical integration integrals
class ANumElement: public AElement{
public:
	~ANumElement() = default;
	ANumElement(const std::vector<int>& nbases, const std::vector<Point>& coo, const IQuadratureRule* quad);
	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;
	std::vector<double> load() const override;
private:
	virtual std::array<double, 9> jacobi_matrix(Point p) const = 0;
	virtual std::vector<double> bases(Point p) const = 0;
	virtual std::vector<Point> grad_bases(Point p) const = 0;
	const IQuadratureRule* _quad;
};

#endif
