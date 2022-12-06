#include "debug.hpp"
#include "prob/poisson_solver.hpp"

PoissonSolver::PoissonSolver(std::shared_ptr<ASpatialApproximator> appr): _approximator(appr) {
	_slae_solver.reset(new AmgcMatrixSolver());
}

void PoissonSolver::set_bc_dirichlet(int btype, double value){
	set_bc_dirichlet(btype, [value](Point){ return value; });
}

void PoissonSolver::set_bc_neumann(int btype, double value){
	set_bc_neumann(btype, [value](Point){ return value; });
}

void PoissonSolver::set_bc_dirichlet(int btype, std::function<double(Point)> value){
	_bc_dirichlet[btype] = value;
}

void PoissonSolver::set_bc_neumann(int btype, std::function<double(Point)> value){
	_bc_neumann[btype] = value;
}

void PoissonSolver::set_bc_robin(int btype, std::function<double(Point)> alpha, std::function<double(Point)> beta){
	_bc_robin_alpha[btype] = alpha;
	_bc_robin_beta[btype] = beta;
}

void PoissonSolver::initialize(){
	_slae_rhs.resize(_approximator->n_bases());
	_rhs_mat = _approximator->mass();
	std::vector<double> diffmat = _approximator->stiff();

	// linear term
	if (_linear_terms.size() > 0){
		// 1. build vector from linear terms
		std::vector<double> lin_vec(_approximator->n_bases(), 0.0);
		for (auto fun: _linear_terms){
			std::vector<double> lin_vec_i = _approximator->approximate(fun);
			for (int i=0; i<(int)lin_vec.size(); ++i){
				lin_vec[i] += lin_vec_i[i];
			}
		}
		// 2. build mass matrix with respect to the linear term
		std::vector<double> mass_alpha = _approximator->mass(lin_vec);
		// 3. add to lhs
		for (size_t i=0; i<diffmat.size(); ++i) diffmat[i] += mass_alpha[i];
	}

	// robin alpha
	for (auto& it: _bc_robin_alpha){
		_approximator->apply_bc_robin_to_stiff_lhs(it.first, it.second, diffmat);
	}

	// dirichlet
	for (auto& it: _bc_dirichlet){
		_approximator->apply_bc_dirichlet_lhs(it.first, diffmat);
	}

	_slae_solver->set_matrix(_approximator->stencil(), diffmat);


	std::cout << "Poisson solver initialized" << std::endl;
}

void PoissonSolver::solve(const std::vector<double>& rhs, std::vector<double>& u){
	// right hand side
	_approximator->stencil().matvec(_rhs_mat, rhs, _slae_rhs);

	// neumann values
	for (auto& it: _bc_neumann){
		_approximator->apply_bc_neumann_to_stiff(it.first, it.second, _slae_rhs);
	}

	// robin beta
	for (auto& it: _bc_robin_beta){
		_approximator->apply_bc_robin_to_stiff_rhs(it.first, it.second, _slae_rhs);
	}
	
	// dirichlet values
	for (auto& it: _bc_dirichlet){
		_approximator->apply_bc_dirichlet_rhs(it.first, it.second, _slae_rhs);
	}

	// solve
	_slae_solver->solve(_slae_rhs, u);

	std::cout << "Poisson solution obtained" << std::endl;
}

// add +eps*u to the equation
void PoissonSolver::add_linear_term(std::function<double(Point)> alpha){
	_linear_terms.push_back(alpha);
}
