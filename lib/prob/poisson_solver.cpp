#include "prob/poisson_solver.hpp"

PoissonSolver::PoissonSolver(std::shared_ptr<ASpatialApproximator> appr): _approximator(appr) {
	_slae_solver.reset(new AmgcMatrixSolver());
}

void PoissonSolver::set_bc_dirichlet(int btype, double value){
	set_bc_dirichlet(btype, [value](Point){ return value; });
}

void PoissonSolver::set_bc_dirichlet(int btype, std::function<double(Point)> value){
	_bc_dirichlet[btype] = value;
}

void PoissonSolver::initialize(){
	_slae_rhs.resize(_approximator->n_bases());
	_rhs_mat = _approximator->mass();
	std::vector<double> diffmat = _approximator->stiff();

	for (auto& it: _bc_dirichlet){
		_approximator->apply_bc_dirichlet_lhs(it.first, diffmat);
	}

	_slae_solver->set_matrix(_approximator->stencil(), diffmat);

	std::cout << "Poisson solver initialized" << std::endl;
}

void PoissonSolver::solve(const std::vector<double>& rhs, std::vector<double>& u){
	// right hand side
	_approximator->stencil().matvec(_rhs_mat, rhs, _slae_rhs);

	// dirichlet values
	for (auto& it: _bc_dirichlet){
		_approximator->apply_bc_dirichlet_rhs(it.first, it.second, _slae_rhs);
	}

	// solve
	_slae_solver->solve(_slae_rhs, u);

	std::cout << "Poisson solution obtained" << std::endl;
}
