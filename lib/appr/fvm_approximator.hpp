#ifndef FVM_APPROXIMATION_HPP
#define FVM_APPROXIMATION_HPP

#include "common.hpp"
#include "grid/regular_grid.hpp"
#include "appr/spatial_approximator.hpp"

class FvmApproximator: public ASpatialApproximator{
public:
	static std::shared_ptr<FvmApproximator> build(std::shared_ptr<AGrid> grid);

	// =============== approximation dimension
	int n_bases() const override;

	// =============== approximate function
	std::vector<double> approximate(std::function<double(Point)> func) const override;

	std::vector<double> mass() const override;
	std::vector<double> stiff() const override;

	void vtk_save_scalars(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const override;
private:
	FvmApproximator(std::shared_ptr<AGrid> grid);
	std::shared_ptr<AGrid> _grid;

	CsrStencil _build_stencil() const override;
	std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const override;
	std::vector<double> _build_load_vector() const override;
};
#endif
