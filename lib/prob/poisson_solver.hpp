#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "common.hpp"
#include "slae/matrix_solver.hpp"
#include "appr/spatial_approximator.hpp"

class PoissonSolver{
public:
	PoissonSolver(std::shared_ptr<ASpatialApproximator> appr);

	void set_bc_dirichlet(int btype, double value);
	void set_bc_dirichlet(int btype, std::function<double(Point)> value);

	void initialize();

	void solve(const std::function<double(Point)>& f, std::vector<double>& u);
	void solve(const std::vector<double>& f, std::vector<double>& u);
private:
	std::shared_ptr<ASpatialApproximator> _approximator;

	std::vector<double> _rhs_mat;
	std::vector<double> _slae_rhs;

	std::unique_ptr<AmgcMatrixSolver> _slae_solver;

	std::map<int, std::function<double(Point)>> _bc_dirichlet;
};

#endif
