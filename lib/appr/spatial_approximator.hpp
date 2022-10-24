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

	// =============== integration
	const std::vector<double>& load_vector() const;
	double integrate(const std::vector<double>& v) const;
	double domain_volume() const;

	// =============== norms
	virtual double norm_max(const std::vector<double>& v) const;
	virtual double norm_2(const std::vector<double>& v) const;
	double rms(const std::vector<double>& v) const;
	double rms(const std::vector<double>& a, const std::vector<double>& b) const;
	double rms(const std::vector<double>& a, std::function<double(Point)> b) const;

	// =============== boundary conditions
	// get boundary basis indices with respective point coordinates
	const std::vector<std::pair<int, Point>>& boundary_bases(int ibnd) const;

	// forced dirichlet conditions
	void apply_bc_dirichlet(
		int ibnd,
		std::function<double(Point)> val_func,
		std::vector<double>& lhs, std::vector<double>& rhs) const;
	void apply_bc_dirichlet_lhs(int ibnd, std::vector<double>& lhs) const;
	void apply_bc_dirichlet_rhs(int ibnd, std::function<double(Point)> val_func, std::vector<double>& rhs) const;

	// stiff matrix conditions
	// du/dn = -q
	virtual void apply_bc_neumann_to_stiff(int ibnd, std::function<double(Point)> q_func, std::vector<double>& rhs) const;
	// du/dn = -alpha*u + beta
	virtual void apply_bc_robin_to_stiff(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::function<double(Point)> beta_func,
		std::vector<double>& stiff, std::vector<double>& rhs) const;

	// =============== savers
	// save single scalar
	void vtk_save_scalar(std::string filepath, const std::vector<double>& scalar, std::string datacap) const;
	// save multiple scalars
	void vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const;
private:
	struct Cache{
		CsrStencil stencil;
		std::vector<double> load_vector;
		double volume = 0;
		std::map<int, std::vector<std::pair<int, Point>>> boundary_bases;
	};
	mutable Cache _cache;

	virtual CsrStencil _build_stencil() const = 0;
	virtual std::vector<double> _build_load_vector() const = 0;
	virtual std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const = 0;

	virtual void _vtk_save_scalar(std::string filepath, std::map<std::string, const std::vector<double>*> scalars) const = 0;
};

#endif
