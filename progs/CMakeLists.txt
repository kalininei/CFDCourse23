# ============= sandbox
add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox PUBLIC ${CFDLIB_TARGET})

# ============= poisson
add_executable(poisson poisson.cpp)
target_link_libraries(poisson PUBLIC ${CFDLIB_TARGET})
