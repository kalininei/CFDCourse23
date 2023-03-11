#ifndef LINEAR_FEM_APPROXIMATION_HPP
#define LINEAR_FEM_APPROXIMATION_HPP

#include "appr/spatial_approximator.hpp"
#include "grid/agrid.hpp"
#include "appr/fe/aelement.hpp"

class LinearFemApproximator: public ASpatialApproximator{
public:
	static std::shared_ptr<LinearFemApproximator> build(std::shared_ptr<AGrid> grid);

	int n_bases() const override;
	std::vector<double> approximate(std::function<double(Point)> func) const override;

	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;

	void apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const override;
protected:
	LinearFemApproximator(std::shared_ptr<AGrid> grid);

	CsrStencil _build_stencil() const override;
	std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const override;
	void vtk_save_scalars(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const override;
	std::vector<double> _build_load_vector() const override;

	int n_elements() const;
	void add_local_matrix(double coef, const std::vector<double>& lmat, const AElement* elem, std::vector<double>& gmat) const;
	void add_local_vector(double coef, const std::vector<double>& lvec, const AElement* elem, std::vector<double>& gvec) const;
	AElement* build_element(int icell);
	AElement* build_boundary_element(int iface);

	std::vector<std::shared_ptr<AElement>> _elements;
	// iface -> element
	std::map<int, std::shared_ptr<AElement>> _boundary_elements;

	const std::shared_ptr<AGrid> _grid;
};


#endif