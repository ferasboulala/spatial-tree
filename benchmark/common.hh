#pragma once

#include <random>
#include <vector>

template <typename CoordT>
std::vector<std::pair<CoordT, CoordT>> generate_points(CoordT BEG, CoordT END, uint64_t test_size) {
    using DistributionType = std::conditional<std::is_floating_point_v<CoordT>,
                                              std::uniform_real_distribution<CoordT>,
                                              std::uniform_int_distribution<CoordT>>::type;
    DistributionType                       distribution(BEG, END);
    std::random_device                     device;
    std::vector<std::pair<CoordT, CoordT>> points;
    for (uint64_t i = 0; i < test_size; ++i) {
        CoordT x = distribution(device);
        CoordT y = distribution(device);
        points.emplace_back(x, y);
    }

    return points;
}