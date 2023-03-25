#ifndef QUADRATIC_FEM_APPROXIMATION_HPP
#define QUADRATIC_FEM_APPROXIMATION_HPP

#include "appr/anodal_fem_approximator.hpp"

class QuadraticFemApproximator: public ANodalFemApproximator{
public:
	static std::shared_ptr<QuadraticFemApproximator> build(std::shared_ptr<AGrid> grid);
protected:
	QuadraticFemApproximator(std::shared_ptr<AGrid> grid);
};

class QuadraticFemApproximator_1D: public QuadraticFemApproximator{
	friend class QuadraticFemApproximator;
	int n_bases() const override;
protected:
	QuadraticFemApproximator_1D(std::shared_ptr<AGrid> grid);

	AElement* build_element(int icell) override;
	AElement* build_boundary_element(int iface) override;

	Point node(int inode) const override;
};

#endif
