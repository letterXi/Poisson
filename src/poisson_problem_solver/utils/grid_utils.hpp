#include <vector>

inline std::vector<size_t> create_block_mask(size_t N, size_t rows, size_t cols) {
    std::vector<size_t> mask(N * N);
    for (size_t j = 0; j < N; ++j) {
        for (size_t i = 0; i < N; ++i) {
            size_t r = (j * rows) / N;
            size_t c = (i * cols) / N;
            if (r >= rows)
                r = rows - 1;
            if (c >= cols)
                c = cols - 1;
            mask[i + j * N] = r * cols + c;
        }
    }
    return mask;
}
