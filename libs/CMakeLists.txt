set(SRC_amgsollib
    src/CsrMatrix.cpp
    src/AmgclSolver.cpp
)

add_library(
    amgsollib
    STATIC
    ${SRC_amgsollib}
)

add_library(
    poissonlib
    STATIC
    src/poisson.cpp
)

target_include_directories(amgsollib PUBLIC include)
target_include_directories(poissonlib PUBLIC include)
