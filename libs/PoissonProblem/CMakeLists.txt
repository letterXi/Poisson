set(SRC
    matrix/src/CsrMatrix.cpp
    matrix_solver/src/AmgclSolver.cpp
    Poisson/src/Poisson.cpp
    shwartz/src/shwartz.cpp
) 

set(HEADERS
    matrix/include/CsrMatrix.hpp
    matrix_solver/include/AmgclSolver.hpp
    Poisson/include/Poisson.hpp
    shwartz/include/shwartz.hpp
)

add_library(PoissonProblemLib STATIC ${SRC} ${HEADERS})

target_include_directories(PoissonProblemLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/matrix/include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/matrix_solver/include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/Poisson/include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/shwartz/include)

                                                #add_subdirectory(test)
