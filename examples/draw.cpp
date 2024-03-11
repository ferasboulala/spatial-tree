#include <raylib.h>

#include <algorithm>

#include "spatial-tree.h"

int main() {
    static constexpr uint64_t WindowWidth = 1200;
    static constexpr uint64_t WindowHeight = 1200;

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
            DrawCircle(x, y, 2, GREEN);
        }
        quadtree.walk([&](auto bbox, bool terminal) {
            if (terminal)
                return;

            auto [x, y] = bbox.origin();
            DrawLine(x, bbox.starts[1], x, bbox.stops[1], GRAY);
            DrawLine(bbox.starts[0], y, bbox.stops[0], y, GRAY);
        });
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
