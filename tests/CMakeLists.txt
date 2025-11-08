add_executable(
    test
    csr_matrix_test.cpp
    solver_test.cpp
    poisson_test.cpp
)

target_include_directories(test PRIVATE ${CMAKE_SOURCE_DIR}/)

target_link_libraries(
    test
    PRIVATE amgsollib
            Catch2::Catch2WithMain
            poissonlib
)
