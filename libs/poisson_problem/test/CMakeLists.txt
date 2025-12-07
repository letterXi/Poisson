add_executable(
    PoissonProblem_test
    matrix_test.cpp
    matrix_solver_test.cpp
    Poisson_test.cpp
)

target_link_libraries(
    PoissonProblem_test
    PRIVATE PoissonProblemLib
            Catch2::Catch2WithMain
)
