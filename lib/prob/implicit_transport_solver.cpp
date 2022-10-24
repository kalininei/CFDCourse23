#include "prob/implicit_transport_solver.hpp"

ImplicitTransportSolver::ImplicitTransportSolver(std::shared_ptr<ASpatialApproximator> appr):
		ANonstationaryProblem(){
	_appr = appr;
}

void ImplicitTransportSolver::set_initial_conditions(std::function<double(Point)> func){
	_u = _appr->approximate(func);
	_uold = _u;
}

void ImplicitTransportSolver::set_velocity(std::function<Vector(Point)> func){
	_vx = _appr->approximate( [&](Point p)->double { return func(p).x; });
	_vy = _appr->approximate( [&](Point p)->double { return func(p).y; });
	_vz = _appr->approximate( [&](Point p)->double { return func(p).z; });

	_need_reassemble = true;
}

const std::vector<double>& ImplicitTransportSolver::u() const{
	return _u;
}

void ImplicitTransportSolver::set_tau(double tau){
	_tau = tau;
	_need_reassemble = true;
}

double ImplicitTransportSolver::_compute_tau(){
	return _tau;
}

void ImplicitTransportSolver::_reassemble(){
	// transport matrix
	_transport_mat = _appr->transport_upwind(_vx, _vy, _vz);
	// lhs matrix = I + tau*K
	std::vector<double> lhs(_mass_mat.size());
	for (size_t i=0; i<lhs.size(); ++i){
		lhs[i] = _mass_mat[i] + _tau * _transport_mat[i];
	}
	// matrix solver
	_mat_solver->set_matrix(_appr->stencil(), lhs);

	_need_reassemble = false;
}

void ImplicitTransportSolver::_solve_next_step(double tau){
	if (_need_reassemble) _reassemble();
	const CsrStencil& s = _appr->stencil();
	std::swap(_u, _uold);
	std::vector<double> rhs(_u.size());
	s.matvec(_mass_mat, _uold, rhs);
	_mat_solver->solve(rhs, _u);
}

void ImplicitTransportSolver::_initialize(){
	_mat_solver.reset(new AmgcMatrixSolver());
	_mass_mat = _appr->mass();
}
