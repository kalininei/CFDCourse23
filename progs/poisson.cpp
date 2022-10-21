#include <iostream>
#include "common.hpp"
#include "prog_common.hpp"
#include "prob/poisson_solver.hpp"
#include "grid/regular_grid.hpp"
#include "appr/fdm_approximator.hpp"

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
	appr->vtk_save_scalar(from_output_path("poisson_fdm.vtk"), x, "data");
}

int main(){
	try{
		fdm2();
		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
