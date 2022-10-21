#ifndef FDM_APPROXIMATION_HPP
#define FDM_APPROXIMATION_HPP

#include "common.hpp"
#include "grid/regular_grid.hpp"
#include "appr/spatial_approximator.hpp"

class FdmApproximator: public ASpatialApproximator{
public:
	static std::shared_ptr<FdmApproximator> build(std::shared_ptr<RegularGrid> grid);
	
	int n_bases() const override;
	std::vector<double> approximate(std::function<double(Point)> func) const override;

	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;
	void apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff) const override;
	void apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> func, std::vector<double>& vec) const override;

	void vtk_save_scalar(std::string filepath, const std::vector<double>& scalar, std::string datacap) const override;
private:
	FdmApproximator(std::shared_ptr<RegularGrid> grid);

	CsrStencil _build_stencil() const override;

	std::shared_ptr<RegularGrid> _grid;
};
#endif
