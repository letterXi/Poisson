set(SRC
	src/CsrMatrix.cpp
	src/AmgclSolver.cpp
)

add_library(amgsollib STATIC ${SRC})
target_include_directories(amgsollib PUBLIC include)


