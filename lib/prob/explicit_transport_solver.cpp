#include "prob/explicit_transport_solver.hpp"

ExplicitTransportSolver::ExplicitTransportSolver(std::shared_ptr<ASpatialApproximator> appr):
		ANonstationaryProblem(){
	_appr = appr;
}

void ExplicitTransportSolver::set_initial_conditions(std::function<double(Point)> func){
	_u = _appr->approximate(func);
	_uold = _u;
}

void ExplicitTransportSolver::set_velocity(std::function<Vector(Point)> func){
	std::vector<double> _vx = _appr->approximate( [&](Point p)->double { return func(p).x; });
	std::vector<double> _vy = _appr->approximate( [&](Point p)->double { return func(p).y; });
	std::vector<double> _vz = _appr->approximate( [&](Point p)->double { return func(p).z; });

	_transport_mat = _appr->transport_upwind(_vx, _vy, _vz);
}

const std::vector<double>& ExplicitTransportSolver::u() const{
	return _u;
}

void ExplicitTransportSolver::set_tau(double tau){
	_tau = tau;
}

double ExplicitTransportSolver::_compute_tau(){
	return _tau;
}

void ExplicitTransportSolver::_solve_next_step(double tau){
	const CsrStencil& s = _appr->stencil();
	std::swap(_u, _uold);

	// u = -tau * K*u_old + I*u_old
	for (int irow=0; irow<_appr->n_bases(); ++irow){
		_u[irow] = (-tau * s.matvec_irow(irow, _transport_mat, _uold)
		            + s.matvec_irow(irow, _mass_mat, _uold));
	}
}

void ExplicitTransportSolver::_initialize(){
	_mass_mat = _appr->mass();
}
