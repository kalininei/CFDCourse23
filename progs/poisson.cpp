#include <iostream>
#include "common.hpp"
#include "prog_common.hpp"
#include "prob/poisson_solver.hpp"
#include "grid/regular_grid.hpp"
#include "appr/fdm_approximator.hpp"

constexpr double PI2 = 8*std::atan(1.0);

double exact_solution(double x){
	return sin(PI2*2*x) + 0.4*cos(PI2*7*x);
}

double exact_solution_d1(double x){
	return 2*PI2*cos(PI2*2*x) - 0.4*7*PI2*sin(PI2*7*x);
}

double exact_solution_d2(double x){
	return -2*PI2*2*PI2*sin(PI2*2*x) - 0.4*7*7*PI2*PI2*cos(PI2*7*x);
}

double rhs_fun(Point p){
	return -exact_solution_d2(p.x);
}

void fdm_poisson(){
	int n_cells = 100;
	// grid
	std::shared_ptr<RegularGrid> grid(new RegularGrid(1, n_cells + 1));
	grid->define_boundary(1, DirectionCode::X_MINUS);
	grid->define_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, exact_solution(0));

	// 1. Dirichlet
	//slv.set_bc_dirichlet(2, exact_solution(1));
	
	// 2. Neumann
	//auto q = [](Point p)->double { return -exact_solution_d1(p.x); }
	//slv.set_bc_neumann(2, q);

	// 3. Robin
	auto alpha = [](Point)->double { return 1; };
	auto beta = [](Point p)->double { return exact_solution(p.x) + exact_solution_d1(p.x); };
	slv.set_bc_robin(2, alpha, beta);

	// rhs
	std::vector<double> rhs = appr->approximate(rhs_fun);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fdm1.vtk"), x, "data");

	// calculate norm
	std::vector<double> appr_exact = appr->approximate([](Point p){ return exact_solution(p.x); });
	std::vector<double> diff(x.size());
	for (size_t i=0; i<x.size(); ++i) {
		diff[i] = x[i] - appr_exact[i];
	}

	double error = appr->norm_2(diff);
	std::cout << n_cells << " " << error << std::endl;
}

void fdm1(){
	// grid
	std::shared_ptr<RegularGrid> grid(new RegularGrid(10, 100));
	grid->define_boundary(1, DirectionCode::X_MINUS);
	grid->define_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(2, 1);

	// rhs
	std::vector<double> rhs = appr->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fdm1.vtk"), x, "data");
}

void fdm2(){
	// grid
	std::shared_ptr<RegularGrid> grid(new RegularGrid(10, 100, 1, 10, 0, 0));
	grid->define_boundary(1, DirectionCode::X_MINUS);
	grid->define_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(2, 1);

	// rhs
	std::vector<double> rhs = appr->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fdm2.vtk"), x, "data");
}

void fdm3(){
	// grid
	std::shared_ptr<RegularGrid> grid(new RegularGrid(1, 30, 1, 30, 1, 30));
	grid->define_boundary(1, DirectionCode::X_MINUS, [](Point p)->bool{
		if (p.y < 0.3) return false;
		if (p.y > 0.7) return false;
		if (p.z < 0.3) return false;
		if (p.z > 0.7) return false;
		return true;
	});
	grid->define_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(2, 1);

	// rhs
	std::vector<double> rhs = appr->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fdm3.vtk"), x, "data");
}

int main(){
	try{
		fdm_poisson();
		//fdm1();
		//fdm2();
		//fdm3();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
