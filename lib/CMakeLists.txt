set(SRC
	prob/poisson_solver.cpp
	slae/csrmat.cpp
	slae/matrix_solver.cpp
	grid/regular_grid.cpp
	grid/regular_grid_boundary.cpp
	appr/spatial_approximator.cpp
	appr/fdm_approximator.cpp
)

set(HEADERS 
	common.hpp
	geom.hpp
	prob/poisson_solver.hpp
	grid/regular_grid.hpp
	grid/regular_grid_boundary.hpp
	appr/spatial_approximator.hpp
	appr/fdm_approximator.hpp
)

add_library(${CFDLIB_TARGET} STATIC ${SRC} ${HEADERS})
