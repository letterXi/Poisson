add_executable(threadbenchmark
	thread_benchmark.cpp
	${CMAKE_SOURCE_DIR}/create_problem/poisson.cpp
)
add_executable(
    app
    main.cpp
    ${CMAKE_SOURCE_DIR}/create_problem/poisson.cpp
)

target_include_directories(app PRIVATE ${CMAKE_SOURCE_DIR}/create_problem/)
target_include_directories(threadbenchmark PRIVATE ${CMAKE_SOURCE_DIR}/create_problem/)
target_link_libraries(app amgsollib) 
