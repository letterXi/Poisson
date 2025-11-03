cmake_minimum_required(VERSION 3.28.3)

project(template)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Wconversion
    -Wunused-const-variable
    -Wlogical-op
    -Wduplicated-cond
    -Wredundant-decls
    -Wnon-virtual-dtor
    -Woverloaded-virtual
    -Wsuggest-override
    -Wzero-as-null-pointer-constant
    -Wsign-promo
    -Wfloat-equal
    -Wcast-qual
    -Wcast-align
    -Wmissing-field-initializers
    -Werror
)

find_package(
    Boost
    REQUIRED
)
find_package(
    Catch2
    REQUIRED
)

FetchContent_Declare(
    amgcl
    GIT_REPOSITORY https://github.com/ddemidov/amgcl
    GIT_TAG 1.4.2
)
FetchContent_GetProperties(amgcl)
if(NOT
   amgcl_POPULATED
)
    FetchContent_Populate(amgcl)
endif()

include_directories(
    SYSTEM
    ${amgcl_SOURCE_DIR}
)

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(libs)
add_subdirectory(tests)
add_subdirectory(app)

