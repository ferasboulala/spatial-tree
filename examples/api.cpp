#include <iostream>

#include "spatial-tree.h"

int main() {
    st::spatial_set<double, 2> my_quad_tree;
    my_quad_tree.emplace({0.0, 0.0});
    my_quad_tree.emplace({1.0, 1.0});

    auto it = my_quad_tree.find({0.0, 0.0});
    if (it == my_quad_tree.end()) {
        std::cerr << "Could not find a point" << std::endl;
        return 1;
    }

    it = my_quad_tree.find({0.5, 0.5});
    if (it != my_quad_tree.end()) {
        std::cerr << "Found a point that should not be there" << std::endl;
        return 2;
    }

    int counter = 0;
    my_quad_tree.find({-5.0, -5.0, 5.0, 5.0}, [&](auto it) {
        ++counter;
        auto [x, y] = *it;
        std::cout << x << ", " << y << std::endl;
    });
    if (counter != 2) {
        std::cerr << "Could not find all points that were added" << std::endl;
        return 3;
    }

    const auto nearest_points = my_quad_tree.nearest({0.1, 0.1}, 1);
    const auto& coordinates = *nearest_points.front();
    if (coordinates[0] != 0 || coordinates[1] != 0) {
        std::cerr << "Should be <0, 0>" << std::endl;
        return 4;
    }

    if (!my_quad_tree.erase({0.0, 0.0})) {
        std::cerr << "Could not erase <0, 0>" << std::endl;
        return 5;
    }

    std::cout << "Success!" << std::endl;

    return 0;
}
