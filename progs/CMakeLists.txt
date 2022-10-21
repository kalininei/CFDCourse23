# ============= sandbox
add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox PUBLIC ${CFDLIB_TARGET})

# ============= poisson
add_executable(poisson poisson.cpp)
target_link_libraries(poisson PUBLIC ${CFDLIB_TARGET})

# ============= transport
add_executable(transport transport.cpp)
target_link_libraries(transport PUBLIC ${CFDLIB_TARGET})
