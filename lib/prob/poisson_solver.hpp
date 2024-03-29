#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "common.hpp"
#include "slae/matrix_solver.hpp"
#include "appr/spatial_approximator.hpp"

class PoissonSolver{
public:
	PoissonSolver(std::shared_ptr<ASpatialApproximator> appr);

	// add +alpha(x)*u to the equation
	void add_linear_term(std::function<double(Point)> alpha);

	void set_bc_dirichlet(int btype, double value);
	void set_bc_neumann(int btype, double value);

	void set_bc_dirichlet(int btype, std::function<double(Point)> value);
	void set_bc_neumann(int btype, std::function<double(Point)> value);
	void set_bc_robin(int btype, std::function<double(Point)> alpha, std::function<double(Point)> beta);
	void set_source(Point coo, double flowrate);

	void initialize();

	void solve(const std::function<double(Point)>& f, std::vector<double>& u);
	void solve(const std::vector<double>& f, std::vector<double>& u);
private:
	std::shared_ptr<ASpatialApproximator> _approximator;

	std::vector<double> _rhs_mat;
	std::vector<double> _slae_rhs;

	std::unique_ptr<AmgcMatrixSolver> _slae_solver;

	std::map<int, std::function<double(Point)>> _bc_dirichlet;
	std::map<int, std::function<double(Point)>> _bc_neumann;
	std::map<int, std::function<double(Point)>> _bc_robin_alpha;
	std::map<int, std::function<double(Point)>> _bc_robin_beta;
	std::vector<std::pair<Point, double>> _bc_point_sources;

	std::vector<std::function<double(Point)>> _linear_terms;
};

#endif
