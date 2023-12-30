#include <limits>
#include <opencv2/opencv.hpp>

#include "spatial-tree.h"

static constexpr uint64_t                       CANVAS_HEIGHT = 1500;
static constexpr uint64_t                       CANVAS_WIDTH = 2500;
static constexpr int                            POINT_RADIUS = 10;
static cv::Mat                                  canvas;
static constexpr uint8_t                        MAXIMUM_NODE_SIZE = 1;
static st::spatial_set<int, MAXIMUM_NODE_SIZE> *points;
static const char                              *window_name = "canvas";

#include <iostream>
void mouse_callback(int event, int x, int y, int, void *) {
    if (event != cv::EVENT_LBUTTONDOWN) return;

    auto                 res = points->nearest(y, x);
    static constexpr int min_distance_squared = POINT_RADIUS * POINT_RADIUS;
    int                  distance_squared = std::numeric_limits<int>::max();
    if (!res.empty()) {
        auto [y_, x_] = *res.front();
        distance_squared = st::internal::euclidean_distance_squared(x, x_, y, y_);
    }

    if (distance_squared <= min_distance_squared) {
        for (auto it : res) {
            auto [y_, x_] = *it;
            points->erase(y_, x_);
        }
    } else {
        points->emplace(y, x);
    }

    canvas = cv::Mat::zeros(CANVAS_HEIGHT, CANVAS_WIDTH, CV_8UC3);
    points->walk([&](auto boundaries, bool is_leaf) {
        if (is_leaf) return;
        auto [top_x, top_y, bottom_x, bottom_y] = boundaries;
        auto mid_x = (top_x + bottom_x) / 2;
        auto mid_y = (top_y + bottom_y) / 2;

        cv::line(canvas, {bottom_y, mid_x}, {top_y, mid_x}, {128, 128, 128}, 2);
        cv::line(canvas, {mid_y, top_x}, {mid_y, bottom_x}, {128, 128, 128}, 2);
    });
    for (auto [y, x] : *points) {
        cv::circle(canvas, {x, y}, POINT_RADIUS, {0, 255, 0}, cv::FILLED);
    }
    cv::imshow(window_name, canvas);
}

int main() {
    cv::namedWindow(window_name);
    cv::setMouseCallback(window_name, mouse_callback);

    canvas = cv::Mat::zeros(CANVAS_HEIGHT, CANVAS_WIDTH, CV_8UC3);
    points = new st::spatial_set<int, MAXIMUM_NODE_SIZE>({0, CANVAS_WIDTH, CANVAS_HEIGHT, 0});
    do {
        cv::imshow(window_name, canvas);
    } while (cv::waitKey(0) != 113);

    return 0;
}