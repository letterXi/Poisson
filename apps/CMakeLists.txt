add_executable(app main.cpp)
add_executable(test_overlap overlap.cpp)

target_link_libraries(
    app
    poisson_problem_lib
)

target_link_libraries(
    test_overlap
    poisson_problem_lib
)
