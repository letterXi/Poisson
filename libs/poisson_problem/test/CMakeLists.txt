add_executable(
    poisson_problem_test
    csr_mat_solver_test.cpp
    csr_matrix_test.cpp
)

target_link_libraries(
    poisson_problem_test
    PRIVATE poisson_problem_lib
            Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(poisson_problem_test)
