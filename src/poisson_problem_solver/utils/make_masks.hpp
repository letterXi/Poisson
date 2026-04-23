#pragma once
#include "poisson_problem_solver/grid/grid_regular2d.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

inline std::vector<size_t> stripes_mask(const RegularGrid2D& grid, size_t n, double angle_deg) {
    const size_t n_total = grid.npoints();
    std::vector<size_t> mask(n_total);

    // Угол в радианы
    const double angle_rad = angle_deg * M_PI / 180.0;
    const double cos_a = std::cos(angle_rad);
    const double sin_a = std::sin(angle_rad);

    // Временный вектор для хранения проекций всех точек
    std::vector<double> projections(n_total);
    double p_min = std::numeric_limits<double>::max();
    double p_max = std::numeric_limits<double>::lowest();

    // 1. Проецируем точки на направление, перпендикулярное полосам
    for (size_t k = 0; k < n_total; ++k) {
        const Point& p = grid.get_point(k);
        // Формула проекции: x*cos(a) + y*sin(a)
        double proj = p.x * cos_a + p.y * sin_a;
        projections[k] = proj;

        if (proj < p_min)
            p_min = proj;
        if (proj > p_max)
            p_max = proj;
    }

    // 2. Вычисляем ширину полосы в пространстве проекций
    double range = p_max - p_min;
    // Используем небольшой эпсилон, чтобы крайняя точка попала в n-1, а не в n область
    double step = (range > 0) ? (range + 1e-9) / static_cast<double>(n) : 1.0;

    // 3. Заполняем маску
    for (size_t k = 0; k < n_total; ++k) {
        size_t stripe_idx = static_cast<size_t>((projections[k] - p_min) / step);

        // Гарантируем попадание в диапазон [0, n-1]
        if (stripe_idx >= n)
            stripe_idx = n - 1;

        mask[k] = stripe_idx;
    }

    return mask;
}

inline std::vector<size_t> block_mask(const RegularGrid2D& grid, size_t rows, size_t cols) {
    const size_t nx = grid.nx();
    const size_t ny = grid.ny();
    std::vector<size_t> mask(grid.npoints());

    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            // Рассчитываем индексы строки и столбца в маске подобластей
            size_t r = (j * rows) / ny;
            size_t c = (i * cols) / nx;

            // Защита от выхода за границы из-за особенностей целочисленного деления
            if (r >= rows)
                r = rows - 1;
            if (c >= cols)
                c = cols - 1;

            // Индекс в маске соответствует j * nx + i (построчное хранение в RegularGrid2D)
            mask[i + j * nx] = r * cols + c;
        }
    }
    return mask;
}

inline std::vector<size_t> diagonal_mask(const RegularGrid2D& grid) {
    const size_t nx = grid.nx();
    const size_t ny = grid.ny();
    std::vector<size_t> mask(grid.npoints());

    // Соотношение сторон для нормировки координат
    // Используем double, чтобы сравнение было корректным для неквадратных сеток
    double aspect = static_cast<double>(ny - 1) / static_cast<double>(nx - 1);

    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            // Уравнения диагоналей в нормализованном виде:
            // Главная: y = aspect * x
            // Побочная: y = (ny - 1) - aspect * x

            bool above_main = (static_cast<double>(j) > aspect * static_cast<double>(i));
            bool above_back = (static_cast<double>(j) > static_cast<double>(ny - 1) - aspect * static_cast<double>(i));

            size_t subdomain_id;
            if (!above_main && !above_back) {
                subdomain_id = 0; // Нижний треугольник
            } else if (above_main && !above_back) {
                subdomain_id = 1; // Левый треугольник
            } else if (above_main && above_back) {
                subdomain_id = 2; // Верхний треугольник
            } else {
                subdomain_id = 3; // Правый треугольник
            }

            mask[i + j * nx] = subdomain_id;
        }
    }
    return mask;
}

inline std::vector<size_t> nested_rects_mask(const RegularGrid2D& grid, size_t n_layers) {
    const size_t nx = grid.nx();
    const size_t ny = grid.ny();
    std::vector<size_t> mask(grid.npoints());

    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            // Нормализуем координаты i и j к расстоянию от краев (от 0 до 0.5)
            double dist_x = std::min(static_cast<double>(i), static_cast<double>(nx - 1 - i)) / nx;
            double dist_y = std::min(static_cast<double>(j), static_cast<double>(ny - 1 - j)) / ny;

            // "Глубина" точки — это минимальное расстояние до ближайшей границы
            double depth = std::min(dist_x, dist_y);

            // Определяем индекс слоя (0 - самый внешний, n_layers-1 - самый внутренний)
            // Умножаем на 2, так как максимальная глубина в центре = 0.5
            size_t layer_idx = static_cast<size_t>(depth * 2.0 * n_layers);

            if (layer_idx >= n_layers) layer_idx = n_layers - 1;

            mask[i + j * nx] = layer_idx;
        }
    }
    return mask;
}
