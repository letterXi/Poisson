#include "../../examples/test_functions.hpp"
#include <thread>
#include "../../main_test.hpp"
#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include "poisson_problem_solver/utils/make_masks.hpp"

TEST_CASE("rectangle partition", "[performance]") {
    size_t N = 500;
    std::vector<size_t> mask;
    RegularGrid2D grid(0, 0, 1, 1, N, N);

    mask = block_mask(grid, 1, 4);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks1x4");

    mask = block_mask(grid, 4, 1);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks4x1");

    mask = block_mask(grid, 2, 2);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks2x2");

   mask = block_mask(grid, 3, 3);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks3x3");

    mask = block_mask(grid, 3, 3);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks4x4_small");

    mask = diagonal_mask(grid);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "diagonals2x2");

    mask = block_mask(grid, 1, 2);
   std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks1x2");

    mask = block_mask(grid, 2, 1);
    std::this_thread::sleep_for(std::chrono::seconds(20));
    log_convergence_metrics(N, mask, "blocks2x1");
}
