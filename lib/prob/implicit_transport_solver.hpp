#ifndef IMPLICIT_TRANSPORT_SOLVER
#define IMPLICIT_TRANSPORT_SOLVER

#include "common.hpp"
#include "prob/nonstationary_problem.hpp"
#include "appr/spatial_approximator.hpp"
#include "slae/matrix_solver.hpp"


class ImplicitTransportSolver: public ANonstationaryProblem{
public:
	ImplicitTransportSolver(std::shared_ptr<ASpatialApproximator> appr);

	void set_initial_conditions(std::function<double(Point)> func);
	void set_velocity(std::function<Vector(Point)> func);
	void set_tau(double tau);

	const std::vector<double>& u() const;
protected:
	double _compute_tau() override;
	void _solve_next_step(double tau) override;
	void _initialize() override;
private:
	void _reassemble();

	double _tau;
	bool _need_reassemble=true;
	std::vector<double> _u, _uold;
	std::vector<double> _vx, _vy, _vz;
	std::shared_ptr<ASpatialApproximator> _appr;
	std::shared_ptr<AmgcMatrixSolver> _mat_solver;
	std::vector<double> _transport_mat;
	std::vector<double> _mass_mat;
};


#endif
