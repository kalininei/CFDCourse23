#include "prog_common.hpp"
#include "geom.hpp"
#include "grid/regular_grid.hpp"
#include "appr/fdm_approximator.hpp"
#include "prob/explicit_transport_solver.hpp"
#include "prob/implicit_transport_solver.hpp"

void fdm1_exp(){
	Vector v {1.0, 0.0, 0.0};

	auto initial_solution = [](Point p)->double{
		if (p.x < 1) return 0;
		if (p.x > 2) return 0;
		return 1;
	};
	auto exact_solution = [initial_solution, v](Point p, double time)->double{
		return initial_solution(p - time*v);
	};

	// grid
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(10, 101);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);
	
	// solver
	ExplicitTransportSolver slv(appr);

	// initial conditions
	slv.set_initial_conditions(initial_solution);

	// velocity
	slv.set_velocity([v](Point)->Vector{ return v; });

	// time step
	slv.set_tau(0.05);

	// ==================== Monitors
	// reports to console each 20 iterations
	auto console_monitor = std::make_shared<ConsoleIterReport>(10);
	slv.add_monitor(console_monitor);
	// saves solution and exact solution to vtk with dt = 0.1
	auto data_saver = std::make_shared<VtkFieldTimeSaver>(0.05, from_output_path("transport_fdm_exp"), appr.get());
	data_saver->add_fun("u", std::bind(&ExplicitTransportSolver::u, &slv));
	data_saver->add_fun("u_exact", [&]()->std::vector<double>{
		auto f = std::bind(exact_solution, std::placeholders::_1, slv.last_time());
		return appr->approximate(f);
	});
	slv.add_monitor(data_saver);
	// reports error each iteration
	auto func_saver = std::make_shared<FunctionalSaver>(1, from_output_path("exp_fun.csv"));
	func_saver->add_fun("Error-2", [&]()->double{
		auto f = std::bind(exact_solution, std::placeholders::_1, slv.last_time());
		return appr->rms(slv.u(), f);
	});
	slv.add_monitor(func_saver);

	// solve at t = 5
	slv.solve(5);
}

void fdm1_imp(){
	Vector v {1.0, 0.0, 0.0};

	auto initial_solution = [](Point p)->double{
		if (p.x < 1) return 0;
		if (p.x > 2) return 0;
		return 1;
	};
	auto exact_solution = [initial_solution, v](Point p, double time)->double{
		return initial_solution(p - time*v);
	};

	// grid
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(10, 101);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);
	
	// solver
	ImplicitTransportSolver slv(appr);

	// initial conditions
	slv.set_initial_conditions(initial_solution);

	// velocity
	slv.set_velocity([v](Point)->Vector{ return v; });
	
	// time step
	slv.set_tau(0.05);

	// ==================== Monitors
	// reports to console each 20 iterations
	auto console_monitor = std::make_shared<ConsoleIterReport>(10);
	slv.add_monitor(console_monitor);
	// saves solution and exact solution to vtk with dt = 0.05
	auto data_saver = std::make_shared<VtkFieldTimeSaver>(0.05, from_output_path("transport_fdm_imp"), appr.get());
	data_saver->add_fun("u", std::bind(&ImplicitTransportSolver::u, &slv));
	data_saver->add_fun("u_exact", [&]()->std::vector<double>{
		auto f = std::bind(exact_solution, std::placeholders::_1, slv.last_time());
		return appr->approximate(f);
	});
	slv.add_monitor(data_saver);
	// reports error each iteration
	auto func_saver = std::make_shared<FunctionalSaver>(1, from_output_path("imp_fun.csv"));
	func_saver->add_fun("Error-2", [&]()->double{
		auto f = std::bind(exact_solution, std::placeholders::_1, slv.last_time());
		return appr->rms(slv.u(), f);
	});
	slv.add_monitor(func_saver);

	// solve at t = 5
	slv.solve(5);
}

int main(){
	try{
		fdm1_exp();
		fdm1_imp();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
