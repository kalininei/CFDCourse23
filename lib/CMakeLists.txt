SET (HEADERS
      common.hpp
      debug.hpp
      geom.hpp

      appr/fdm_approximator.hpp
      appr/fvm_approximator.hpp
      appr/output_fields.hpp
      appr/spatial_approximator.hpp
      appr/linear_fem_approximator.hpp
      appr/fe/aelement.hpp
      appr/fe/linear_segment.hpp
      appr/fe/linear_triangle.hpp
      appr/fe/linear_tetrahedron.hpp


      grid/agrid.hpp
      grid/agrid_boundary.hpp
      grid/regular_grid.hpp
      grid/regular_grid_boundary.hpp
      grid/unstructured_grid.hpp

      prob/explicit_transport_solver.hpp
      prob/implicit_transport_solver.hpp
      prob/monitor.hpp
      prob/nonstationary_problem.hpp
      prob/poisson_solver.hpp

      slae/csrmat.hpp
      slae/matrix_solver.hpp
)

SET (SRC
      debug.cpp

      appr/fdm_approximator.cpp
      appr/fvm_approximator.cpp
      appr/spatial_approximator.cpp
      appr/linear_fem_approximator.cpp
      appr/fe/aelement.cpp
      appr/fe/linear_segment.cpp
      appr/fe/linear_triangle.cpp
      appr/fe/linear_tetrahedron.cpp

      grid/agrid.cpp
      grid/agrid_boundary.cpp
      grid/regular_grid.cpp
      grid/regular_grid_boundary.cpp
      grid/unstructured_grid.cpp

      prob/poisson_solver.cpp
      prob/nonstationary_problem.cpp
      prob/monitor.cpp
      prob/implicit_transport_solver.cpp
      prob/explicit_transport_solver.cpp

      slae/matrix_solver.cpp
      slae/csrmat.cpp


)


add_library(${CFDLIB_TARGET} STATIC ${HEADERS} ${SRC})
