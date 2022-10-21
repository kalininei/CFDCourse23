#include <iostream>
#include "common.hpp"
#include "prog_common.hpp"
#include "prob/poisson_solver.hpp"
#include "grid/regular_grid.hpp"
#include "appr/fdm_approximator.hpp"

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
		fdm1();
		fdm2();
		fdm3();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
