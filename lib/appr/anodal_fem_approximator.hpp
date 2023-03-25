#ifndef ANODAL_FEM_APPROXIMATION_HPP
#define ANODAL_FEM_APPROXIMATION_HPP

#include "appr/spatial_approximator.hpp"
#include "grid/agrid.hpp"
#include "appr/fe/aelement.hpp"

class ANodalFemApproximator: public ASpatialApproximator{
public:
	std::vector<double> approximate(std::function<double(Point)> func) const override;

	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;

	void apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const override;
	// du/dn = -alpha*u + beta
	void apply_bc_robin_to_stiff_lhs(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::vector<double>& stiff) const override;
	void apply_bc_robin_to_stiff_rhs(
		int ibnd,
		std::function<double(Point)> beta_func,
		std::vector<double>& rhs) const override;
protected:
	ANodalFemApproximator(std::shared_ptr<AGrid> grid);
	void initialize();

	CsrStencil _build_stencil() const override;
	std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const override;
	void vtk_save_scalars(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const override;
	std::vector<double> _build_load_vector() const override;

	int n_elements() const;
	void add_local_matrix(double coef, const std::vector<double>& lmat, const AElement* elem, std::vector<double>& gmat) const;
	void add_local_vector(double coef, const std::vector<double>& lvec, const AElement* elem, std::vector<double>& gvec) const;

	virtual AElement* build_element(int icell) = 0;
	virtual AElement* build_boundary_element(int iface) = 0;
	virtual Point node(int inode) const = 0;

	std::vector<std::shared_ptr<AElement>> _elements;
	// iface -> element
	std::map<int, std::shared_ptr<AElement>> _boundary_elements;

	const std::shared_ptr<AGrid> _grid;
};


#endif
