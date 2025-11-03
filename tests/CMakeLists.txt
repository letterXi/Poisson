add_executable(test
	csr_matrix_test.cpp
	solver_test.cpp
)
target_link_libraries(test PRIVATE amgsollib Catch2::Catch2WithMain)
