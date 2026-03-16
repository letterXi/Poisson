set(SRC
  main_test.cpp
  examples/two_equal_rectangles.cpp
  examples/four_squares.cpp
  examples/ring_subdomain.cpp
  examples/parallel_jacobi.cpp
  examples/overlap_convergence_study.cpp
  subdomain_test.cpp
  schwarz_solver_test.cpp
  mat_test.cpp
  mat_solver_test.cpp
)

set(HEADERS
  main_test.hpp
  examples/test_functions.hpp
)

add_executable(
  main_test
  ${SRC}
  ${HEADERS}
)

target_link_libraries(main_test PRIVATE Catch2::Catch2WithMain)
target_link_libraries(main_test PRIVATE poisson_problem_solver)
target_include_directories(main_test PRIVATE ${CMAKE_SOURCE_DIR}/src/test)
