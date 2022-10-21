#ifndef EXPLICIT_TRANSPORT_SOLVER
#define EXPLICIT_TRANSPORT_SOLVER

#include "common.hpp"
#include "prob/nonstationary_problem.hpp"
#include "appr/spatial_approximator.hpp"


class ExplicitTransportSolver: public ANonstationaryProblem{
public:
	ExplicitTransportSolver(std::shared_ptr<ASpatialApproximator> appr);

	void set_initial_conditions(std::function<double(Point)> func);
	void set_velocity(std::function<Vector(Point)> func);
	void set_tau(double tau);

	const std::vector<double>& u() const;
protected:
	double _compute_tau() override;
	void _solve_next_step(double tau) override;
	void _initialize() override;

private:
	double _tau;
	std::vector<double> _u, _uold;
	std::shared_ptr<ASpatialApproximator> _appr;

	std::vector<double> _transport_mat;
	std::vector<double> _mass_mat;
};


#endif
