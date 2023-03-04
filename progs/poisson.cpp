#include <iostream>
#include "common.hpp"
#include "prog_common.hpp"
#include "prob/poisson_solver.hpp"
#include "grid/regular_grid.hpp"
#include "grid/unstructured_grid.hpp"
#include "appr/fdm_approximator.hpp"
#include "appr/fvm_approximator.hpp"

const double PI2 = 8*std::atan(1.0);

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
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(n_cells + 1, 1);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<FdmApproximator> appr = FdmApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// add +eps*u to the equation
	auto eps = [](Point)->double { return 1e-6; };
	slv.add_linear_term(eps);

	// bc
	// ==== Left side
	//// Dirichlet
	// slv.set_bc_dirichlet(1, exact_solution(0));

	// ... or Neumann
	auto q1 = [](Point p)->double { return exact_solution_d1(p.x); };
	slv.set_bc_neumann(1, q1);

	// ==== right side
	//// Dirichlet
	//slv.set_bc_dirichlet(2, exact_solution(1));
	
	// ... or Neumann
	auto q = [](Point p)->double { return -exact_solution_d1(p.x); };
	slv.set_bc_neumann(2, q);

	//// ... or Robin
	//auto alpha = [](Point)->double { return 1; };
	//auto beta = [](Point p)->double { return exact_solution(p.x) + exact_solution_d1(p.x); };
	//slv.set_bc_robin(2, alpha, beta);

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

void fvm_poisson(){
	int n_cells = 100;
	// grid
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(n_cells + 1, 1);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<ASpatialApproximator> appr = FvmApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	// ==== Left side
	// Dirichlet
	slv.set_bc_dirichlet(1, exact_solution(0));

	// ==== right side
	// Dirichlet
	//slv.set_bc_dirichlet(2, exact_solution(1));
	
	//// Neumann
	//auto q = [](Point p)->double { return -exact_solution_d1(p.x); };
	//slv.set_bc_neumann(2, q);
	
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
	appr->vtk_save_scalar(from_output_path("poisson_fvm1.vtk"), x, "data");

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
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(100, 10);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

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
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(100, 10, 10, 1, 1, 0);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

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
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(30, 1, 30, 1, 30, 1);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS, [](Point p)->bool{
		if (p.y < 0.3) return false;
		if (p.y > 0.7) return false;
		if (p.z < 0.3) return false;
		if (p.z > 0.7) return false;
		return true;
	});
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

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

void fvm2(){
	std::shared_ptr<UnstructuredGrid> grid = UnstructuredGrid::read_from_vtk(from_input_path("rect.vtk"));
	grid->define_boundary(1, [](Point p)->bool {
		if (p.x < 1e-6 && p.y > 0.2 && p.y < 0.8)
			return true;
		else
			return false;
	});

	// spatial approximator
	std::shared_ptr<ASpatialApproximator> appr = FvmApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	// Dirichlet
	slv.set_bc_dirichlet(1, [](Point p)->double{ return 1; });

	// rhs
	std::vector<double> rhs = appr->approximate(rhs_fun);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fvm2.vtk"), x, "data");

	// calculate norm
	std::vector<double> appr_exact = appr->approximate([](Point p){ return exact_solution(p.x); });
	std::vector<double> diff(x.size());
	for (size_t i=0; i<x.size(); ++i) {
		diff[i] = x[i] - appr_exact[i];
	}

	double error = appr->norm_2(diff);
	std::cout << grid->n_cells() << " " << error << std::endl;
}

void fvm2_neumann(){
	std::shared_ptr<UnstructuredGrid> grid = UnstructuredGrid::read_from_vtk(from_input_path("rect.vtk"));
	grid->define_boundary(1, [](Point p)->bool {
		if (p.x < 1e-6 && p.y > 0.5)
			return true;
		else
			return false;
	});
	grid->define_boundary(2, [](Point p)->bool {
		if (p.x > 1 - 1e-6 && p.y < 0.5)
			return true;
		else
			return false;
	});

	// spatial approximator
	std::shared_ptr<ASpatialApproximator> appr = FvmApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	// Dirichlet
	slv.set_bc_dirichlet(1, [](Point p)->double{ return 1; });
	slv.set_bc_dirichlet(2, [](Point p)->double{ return 0; });

	// point source
	slv.set_source(Point{0.5, 0.5}, 1);
	slv.set_source(Point{0.3, 0.2}, -1);

	// rhs
	std::vector<double> rhs(appr->n_bases(), 0.0);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fvm2.vtk"), x, "data");

	// calculate q
	double q = appr->calculate_dudn(1, x);

	std::cout << grid->n_cells() << " " << q << std::endl;
}

int main(){
	try{
		//fdm_poisson();
		//fdm1();
		//fdm2();
		//fdm3();

		//fvm_poisson();
		//fvm2();
		fvm2_neumann();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
