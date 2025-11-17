add_executable(threadbenchmark thread_benchmark.cpp)
add_executable(app main.cpp)

target_link_libraries(
    app
    PoissonProblemLib
)

target_link_libraries(
    threadbenchmark
    PoissonProblemLib
)
