#include <raylib.h>

#include <algorithm>

#include "spatial-tree.h"

#if 0

static constexpr int     LINE_THICKNESS = 1;
static constexpr uint8_t MAXIMUM_NODE_SIZE = 1;
static const char       *window_name = "canvas";

static constexpr int                               CANVAS_HEIGHT = 1500;
static constexpr int                               POINT_RADIUS = 10;
static constexpr int                               CANVAS_WIDTH = 2500;
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

int main() {
    static constexpr uint64_t WindowWidth = 800;
    static constexpr uint64_t WindowHeight = 600;

    st::spatial_set<int32_t, 2, 1> quadtree({0, 0, WindowWidth, WindowHeight});

    InitWindow(WindowWidth, WindowHeight, "draw");
    SetTargetFPS(60);

    BeginDrawing();
    ClearBackground(RAYWHITE);
    EndDrawing();

    while (!WindowShouldClose()) {
        // TODO: Follow the steps here to speed up drawing:
        // https://www.reddit.com/r/raylib/comments/i6mkh0/only_clear_background_once/
        BeginDrawing();
        ClearBackground(RAYWHITE);
        for (auto [x, y] : quadtree) {
            DrawCircle(x, y, 1, GREEN);
        }
        DrawFPS(10, 10);
        std::string size_str = std::to_string(quadtree.size());
        DrawText(size_str.c_str(), 10, 30, 10, DARKGRAY);

        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            auto mouse_position = GetMousePosition();
            if (mouse_position.x >= 0 && mouse_position.y > 0 && mouse_position.x < WindowWidth &&
                mouse_position.y < WindowHeight) {
                quadtree.emplace({int32_t(mouse_position.x), int32_t(mouse_position.y)});
            }
        } else if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            auto                     mouse_position = GetMousePosition();
            static constexpr int32_t EraserRadius = 20;

            int32_t min_x = std::max<int32_t>(0, mouse_position.x - EraserRadius);
            int32_t max_x = std::min<int32_t>(WindowWidth, mouse_position.x + EraserRadius);
            int32_t min_y = std::max<int32_t>(0, mouse_position.y - EraserRadius);
            int32_t max_y = std::min<int32_t>(WindowHeight, mouse_position.y + EraserRadius);

            std::vector<std::array<int32_t, 2> > points_to_erase;
            quadtree.find({min_x, min_y, max_x, max_y},
                          [&](auto it) { points_to_erase.push_back(it); });
            for (auto point : points_to_erase) {
                quadtree.erase(point);
            }
#define LIGHT_BEIGE \
    CLITERAL(Color) { 200, 176, 131, 255 }
            DrawRectangle(min_x, min_y, 2 * EraserRadius, 2 * EraserRadius, LIGHT_BEIGE);
        }
        EndDrawing();
    }

    CloseWindow();
    return 0;
}

#endif
