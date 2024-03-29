SET (HEADERS
      common.hpp
      debug.hpp
      geom.hpp
      tictoc.hpp
      timedep_writer.hpp

      appr/fdm_approximator.hpp
      appr/fvm_approximator.hpp
      appr/output_fields.hpp
      appr/spatial_approximator.hpp
      appr/anodal_fem_approximator.hpp
      appr/linear_fem_approximator.hpp
      appr/quadratic_fem_approximator.hpp
      appr/cubic_fem_approximator.hpp
      appr/fe/fe_quadrature.hpp
      appr/fe/aelement.hpp
      appr/fe/anum_element.cpp
      appr/fe/linear_segment.hpp
      appr/fe/quadratic_segment.hpp
      appr/fe/cubic_segment.hpp
      appr/fe/linear_triangle.hpp
      appr/fe/linear_tetrahedron.hpp
      appr/fe/bilinear_quadrangle.hpp
      appr/fe/linear_bubble_triangle.hpp
      appr/fe/quadratic_segment.hpp

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
      slae/dense_mat.hpp
)

SET (SRC
      debug.cpp
      tictoc.cpp
      timedep_writer.cpp
      
      appr/fdm_approximator.cpp
      appr/fvm_approximator.cpp
      appr/spatial_approximator.cpp
      appr/anodal_fem_approximator.cpp
      appr/linear_fem_approximator.cpp
      appr/quadratic_fem_approximator.cpp
      appr/cubic_fem_approximator.cpp
      appr/fe/fe_quadrature.cpp
      appr/fe/aelement.cpp
      appr/fe/anum_element.cpp
      appr/fe/linear_segment.cpp
      appr/fe/quadratic_segment.cpp
      appr/fe/cubic_segment.cpp
      appr/fe/linear_triangle.cpp
      appr/fe/linear_tetrahedron.cpp
      appr/fe/bilinear_quadrangle.cpp
      appr/fe/linear_bubble_triangle.cpp
      appr/fe/quadratic_segment.cpp

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
      slae/dense_mat.cpp
)


add_library(${CFDLIB_TARGET} STATIC ${HEADERS} ${SRC})
