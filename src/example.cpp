#include <iostream>

#include "spatial-tree.h"

int main() {
    st::spatial_set<double> my_quad_tree;
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
    my_quad_tree.find(st::bounding_box<double, 2>({-5.0, 5.0, 5.0, -5.0}),
                      [&](auto) { ++counter; });
    if (counter != 2) {
        std::cerr << "Could not find all points that were added" << std::endl;
        return 3;
    }

    auto nearest_points = my_quad_tree.nearest({0.1, 0.1});
    if (nearest_points.size() != 1) {
        std::cerr << "Not the right number of points" << std::endl;
        return 4;
    }

    auto coordinates = *nearest_points.front();
    if (coordinates[0] || coordinates[1]) {
        std::cerr << "Should be <0, 0>" << std::endl;
        return 5;
    }

    return 0;
}