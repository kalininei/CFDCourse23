#ifndef LINEAR_FEM_APPROXIMATION_HPP
#define LINEAR_FEM_APPROXIMATION_HPP

#include "appr/anodal_fem_approximator.hpp"

class LinearFemApproximator: public ANodalFemApproximator{
public:
	static std::shared_ptr<LinearFemApproximator> build(std::shared_ptr<AGrid> grid);
	int n_bases() const override;
protected:
	LinearFemApproximator(std::shared_ptr<AGrid> grid);

	AElement* build_element(int icell) override;
	AElement* build_boundary_element(int iface) override;

	Point node(int inode) const override;
};


#endif
