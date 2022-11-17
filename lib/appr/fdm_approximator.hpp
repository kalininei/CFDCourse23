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
	std::vector<double> transport(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const override;
	std::vector<double> transport_upwind(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const override;
private:
	FdmApproximator(std::shared_ptr<RegularGrid> grid);

	CsrStencil _build_stencil() const override;
	std::vector<double> _build_load_vector() const override;
	std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const override;
	void _vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const override;

	std::shared_ptr<RegularGrid> _grid;
};
#endif
