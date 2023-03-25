#ifndef CUBIC_FEM_APPROXIMATION_HPP
#define CUBIC_FEM_APPROXIMATION_HPP

#include "appr/anodal_fem_approximator.hpp"

class CubicFemApproximator: public ANodalFemApproximator{
public:
	static std::shared_ptr<CubicFemApproximator> build(std::shared_ptr<AGrid> grid);
protected:
	CubicFemApproximator(std::shared_ptr<AGrid> grid);
};

class CubicFemApproximator_1D: public CubicFemApproximator{
	friend class CubicFemApproximator;
	int n_bases() const override;
protected:
	CubicFemApproximator_1D(std::shared_ptr<AGrid> grid);

	AElement* build_element(int icell) override;
	AElement* build_boundary_element(int iface) override;

	Point node(int inode) const override;
};

#endif
