#include <limits>
#include <opencv2/opencv.hpp>

#include "spatial-tree.h"

static constexpr int     LINE_THICKNESS = 1;
static constexpr uint8_t MAXIMUM_NODE_SIZE = 1;
static const char       *window_name = "canvas";

#define DRAW
#ifdef DRAW
static constexpr int                               CANVAS_HEIGHT = 1500;
static constexpr int                               POINT_RADIUS = 10;
static constexpr int                               CANVAS_WIDTH = 2500;
static cv::Mat                                     canvas;
static st::spatial_set<int, 2, MAXIMUM_NODE_SIZE> *points;

void mouse_callback(int event, int x, int y, int, void *) {
    if (event != cv::EVENT_LBUTTONDOWN) return;

    auto                 res = points->nearest({y, x});
    static constexpr int min_distance_squared = POINT_RADIUS * POINT_RADIUS;
    int                  distance_squared = std::numeric_limits<int>::max();
    if (!res.empty()) {
        auto [y_, x_] = res.front();
        distance_squared = st::internal::euclidean_distance_squared(x, x_, y, y_);
    }

    if (distance_squared <= min_distance_squared) {
        for (auto pt : res) {
            points->erase(pt);
        }
    } else {
        points->emplace({y, x});
    }

    canvas = cv::Mat::zeros(CANVAS_HEIGHT, CANVAS_WIDTH, CV_8UC3);
    points->walk([&](auto boundaries, bool is_leaf) {
        if (is_leaf) return;
        auto top_x = boundaries.starts[0];
        auto bottom_y = boundaries.starts[1];
        auto bottom_x = boundaries.stops[0];
        auto top_y = boundaries.stops[1];

        auto mid_x = (top_x + bottom_x) / 2;
        auto mid_y = (top_y + bottom_y) / 2;

        cv::line(canvas, {bottom_y, mid_x}, {top_y, mid_x}, {128, 128, 128}, LINE_THICKNESS);
        cv::line(canvas, {mid_y, top_x}, {mid_y, bottom_x}, {128, 128, 128}, LINE_THICKNESS);
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
    st::spatial_set<int, 2, MAXIMUM_NODE_SIZE> tree({0, 0, CANVAS_HEIGHT, CANVAS_WIDTH});
    points = &tree;
    do {
        cv::imshow(window_name, canvas);
    } while (cv::waitKey(0) != 113);

    return 0;
}

#else

#include <string>
int main(int, char **argv) {
    const std::string input_filename(argv[1]);
    const std::string output_filename(argv[2]);
    auto              map = cv::imread(input_filename, cv::IMREAD_GRAYSCALE);
    cv::threshold(map, map, 128, std::numeric_limits<uint8_t>::max(), cv::THRESH_BINARY);

    st::spatial_set<int, MAXIMUM_NODE_SIZE> points({0, map.cols, map.rows, 0});
    for (int i = 0; i < map.rows; ++i) {
        for (int j = 0; j < map.cols; ++j) {
            if (!map.at<uint8_t>(i, j)) {
                points.emplace(i, j);
            }
        }
    }

    cv::Mat color_map;
    cv::cvtColor(map, color_map, cv::COLOR_GRAY2RGB);

    points.walk([&](auto boundaries, bool is_leaf) {
        if (is_leaf) return;
        auto [top_x, top_y, bottom_x, bottom_y] = boundaries;
        auto mid_x = (top_x + bottom_x) / 2;
        auto mid_y = (top_y + bottom_y) / 2;

        cv::line(color_map, {bottom_y, mid_x}, {top_y, mid_x}, {0, 128, 0}, LINE_THICKNESS);
        cv::line(color_map, {mid_y, top_x}, {mid_y, bottom_x}, {0, 128, 0}, LINE_THICKNESS);
    });

    cv::namedWindow(window_name);
    cv::imshow(window_name, color_map);
    cv::waitKey(0);

    cv::imwrite(output_filename, color_map);

    return 0;
}

#endif