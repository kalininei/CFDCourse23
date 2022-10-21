#ifndef SPATIAL_APPROXIMATOR_HPP
#define SPATIAL_APPROXIMATOR_HPP

#include "common.hpp"
#include "geom.hpp"
#include "slae/csrmat.hpp"

class ASpatialApproximator{
public:
	virtual ~ASpatialApproximator() = default;

	// =============== approximation dimension
	virtual int n_bases() const = 0;

	// =============== approximate function
	// approximates analytic function
	virtual std::vector<double> approximate(std::function<double(Point)> func) const = 0;

	// =============== stiff matrix stencil
	const CsrStencil& stencil() const;

	// =============== approximated differential operators
	// returns matrix for Laplace operator with df/dn = 0 boundry condition
	virtual std::vector<double> stiff() const;
	// returns matrix for Identity operator
	virtual std::vector<double> mass() const;
	// returns matrix for Transport operator
	virtual std::vector<double> transport(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;
	// returns matrix for Upwind Transport operator
	virtual std::vector<double> transport_upwind(std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz) const;

	// =============== boundary conditions
	virtual void apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff_mat) const;
	virtual void apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> dir_values,
	                                             std::vector<double>& stiff_vec) const;

	// =============== savers
	// save single scalar
	void vtk_save_scalar(std::string filepath, const std::vector<double>& scalar, std::string datacap) const;
	// save multiple scalars
	void vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const;
protected:
	virtual CsrStencil _build_stencil() const = 0;
private:
	struct Cache{
		CsrStencil stencil;
	};
	mutable Cache _cache;

	virtual void _vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const = 0;
};

#endif
