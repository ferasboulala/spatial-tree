#include <iostream>
#include <opencv2/opencv.hpp>

#include "spatial.h"

static constexpr uint64_t                       CANVAS_HEIGHT = 1000;
static constexpr uint64_t                       CANVAS_WIDTH = 1000;
static cv::Mat                                  canvas;
static constexpr uint8_t                        MAXIMUM_NODE_SIZE = 1;
static st::spatial_set<int, MAXIMUM_NODE_SIZE> *points;
static const char                              *window_name = "canvas";

void mouse_callback(int event, int x, int y, int, void *) {
    if (event != cv::EVENT_LBUTTONDOWN) return;

    std::cout << "Processing callback" << std::endl;
    if (!points->emplace(x, y).second) {
        points->erase(x, y);
    }

    canvas = cv::Mat::zeros(CANVAS_HEIGHT, CANVAS_WIDTH, CV_8UC3);
    points->walk([&](auto boundaries) {
        auto [top_x, top_y, bottom_x, bottom_y] = boundaries;
        auto mid_x = (top_x + bottom_x) / 2;
        auto mid_y = (top_y + bottom_y) / 2;
        std::cout << mid_y << std::endl;

        cv::line(canvas, {mid_x, bottom_y}, {mid_x, top_y}, {128, 128, 128}, 3);
        cv::line(canvas, {top_x, mid_y}, {bottom_x, mid_y}, {128, 128, 128}, 3);
    });
    for (auto [x, y] : *points) {
        cv::circle(canvas, {x, y}, 5, {0, 255, 0}, cv::FILLED);
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