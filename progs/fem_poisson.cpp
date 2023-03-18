#include <iostream>
#include "prog_common.hpp"
#include "appr/linear_fem_approximator.hpp"
#include "grid/regular_grid.hpp"
#include "grid/unstructured_grid.hpp"
#include "prob/poisson_solver.hpp"

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


void linear1(){
	int n_cells = 10;

	// grid
	std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(n_cells + 1, 1);
	grid->define_reg_boundary(1, DirectionCode::X_MINUS);
	grid->define_reg_boundary(2, DirectionCode::X_PLUS);

	// spatial approximator
	std::shared_ptr<LinearFemApproximator> appr = LinearFemApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, exact_solution(0));
	slv.set_bc_dirichlet(2, exact_solution(1));
	// slv.set_bc_neumann(2, -exact_solution_d1(1));
	
	// rhs
	std::vector<double> rhs = appr->approximate(rhs_fun);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fem1.vtk"), x, "data");

	// calculate norm
	std::vector<double> appr_exact = appr->approximate([](Point p){ return exact_solution(p.x); });
	std::vector<double> diff(x.size());
	for (size_t i=0; i<x.size(); ++i) {
		diff[i] = x[i] - appr_exact[i];
	}

	double error = appr->norm_2(diff);
	std::cout << n_cells << " " << error << std::endl;
};

void linear2(){
	// grid
	std::shared_ptr<UnstructuredGrid> grid = UnstructuredGrid::read_from_vtk(from_input_path("rect.vtk"));
	grid->define_boundary(1, [](Point p)->bool {
		if (p.x < 1e-6)
			return true;
		else
			return false;
	});
	grid->define_boundary(2, [](Point p)->bool {
		if (p.x > 1 - 1e-6)
			return true;
		else
			return false;
	});

	// spatial approximator
	std::shared_ptr<LinearFemApproximator> appr = LinearFemApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, [](const Point& p)->double{ return exact_solution(p.x); });
	slv.set_bc_dirichlet(2, [](const Point& p)->double{ return exact_solution(p.x); });
	//slv.set_bc_neumann(2, [](const Point& p)->double{ return exact_solution_d1(p.x); });
	
	// rhs
	std::vector<double> rhs = appr->approximate(rhs_fun);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fem2.vtk"), x, "data");

	// calculate norm
	std::vector<double> appr_exact = appr->approximate([](Point p){ return exact_solution(p.x); });
	std::vector<double> diff(x.size());
	for (size_t i=0; i<x.size(); ++i) {
		diff[i] = x[i] - appr_exact[i];
	}

	double error = appr->norm_2(diff);
	std::cout << grid->n_cells() << " " << error << std::endl;
};

void linear3(){
	// grid
	std::shared_ptr<UnstructuredGrid> grid = UnstructuredGrid::read_from_vtk(from_input_path("cube.vtk"));
	grid->define_boundary(1, [](Point p)->bool {
		if (p.x < 1e-6)
			return true;
		else
			return false;
	});
	grid->define_boundary(2, [](Point p)->bool {
		if (p.x > 1 - 1e-6)
			return true;
		else
			return false;
	});

	// spatial approximator
	std::shared_ptr<LinearFemApproximator> appr = LinearFemApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, [](const Point& p)->double{ return exact_solution(p.x); });
	slv.set_bc_dirichlet(2, [](const Point& p)->double{ return exact_solution(p.x); });
	
	// rhs
	std::vector<double> rhs = appr->approximate(rhs_fun);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fem3.vtk"), x, "data");

	// calculate norm
	std::vector<double> appr_exact = appr->approximate([](Point p){ return exact_solution(p.x); });
	std::vector<double> diff(x.size());
	for (size_t i=0; i<x.size(); ++i) {
		diff[i] = x[i] - appr_exact[i];
	}

	double error = appr->norm_2(diff);
	std::cout << grid->n_cells() << " " << error << std::endl;
};

void bilinear2(){
	// grid
	//int ndim = 10;
	//std::shared_ptr<ARegularGrid> grid = ARegularGrid::build(ndim + 1, 1, ndim + 1, 1);
	std::shared_ptr<UnstructuredGrid> grid = UnstructuredGrid::read_from_vtk(from_input_path("rect_4_4.vtk"));
	grid->define_boundary(1, [](Point p)->bool {
		if (p.x < 1e-6)
			return true;
		else
			return false;
	});
	grid->define_boundary(2, [](Point p)->bool {
		if (p.x > 1 - 1e-6)
			return true;
		else
			return false;
	});

	// spatial approximator
	std::shared_ptr<LinearFemApproximator> appr = LinearFemApproximator::build(grid);

	// solver: -Laplace(u) = f
	PoissonSolver slv(appr);

	// bc
	slv.set_bc_dirichlet(1, [](const Point& p)->double{ return exact_solution(p.x); });
	slv.set_bc_dirichlet(2, [](const Point& p)->double{ return exact_solution(p.x); });
	
	// rhs
	std::vector<double> rhs = appr->approximate(rhs_fun);

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fem2.vtk"), x, "data");

	// calculate norm
	std::vector<double> appr_exact = appr->approximate([](Point p){ return exact_solution(p.x); });
	std::vector<double> diff(x.size());
	for (size_t i=0; i<x.size(); ++i) {
		diff[i] = x[i] - appr_exact[i];
	}

	double error = appr->norm_2(diff);
	std::cout << grid->n_cells() << " " << error << std::endl;
};

int main(){
	try{
		//linear1();
		linear2();
		//linear3();
		//bilinear2();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
