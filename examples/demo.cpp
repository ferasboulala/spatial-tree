#include <raylib.h>
#include <rlImGui.h>

#include <algorithm>

#include "spatial-tree.h"

class my_spatial_set : public st::spatial_set<int32_t, 2> {
public:
    void draw() const {
        const std::function<void(int32_t, const st::bounding_box<int32_t, 2>&)> walk_recursively =
            [&](auto branch_index, auto boundary) {
                assert(branch_index < this->branches_.size());
                auto& branch = this->branches_[branch_index];

                if (branch.is_terminal()) {
                    return;
                }

                auto [x, y] = boundary.origin();
                DrawLine(x, boundary.starts[1], x, boundary.stops[1], GRAY);
                DrawLine(boundary.starts[0], y, boundary.stops[0], y, GRAY);

                st::internal::unroll_for<BranchingFactor>([&](auto child) {
                    const uint64_t child_index = branch.index_of_first_child + child;
                    const auto     new_boundary = boundary.qrecurse(child);
                    walk_recursively(child_index, new_boundary);
                });
            };

        walk_recursively(0, this->boundary_);
    }
};

int main() {
    static constexpr uint64_t WindowWidth = 800;
    static constexpr uint64_t WindowHeight = 800;

    my_spatial_set quadtree({0, 0, WindowWidth, WindowHeight});

    InitWindow(WindowWidth, WindowHeight, "draw");
    rlImGuiSetup(false);

    BeginDrawing();
    ClearBackground(RAYWHITE);
    EndDrawing();

    while (!WindowShouldClose()) {
        // TODO: Follow the steps here to speed up drawing:
        // https://www.reddit.com/r/raylib/comments/i6mkh0/only_clear_background_once/
        BeginDrawing();
        rlImGuiBegin();
        ClearBackground(RAYWHITE);
        for (auto [x, y] : quadtree) {
            DrawCircle(x, y, 2, GREEN);
        }
        quadtree.draw();
        DrawFPS(10, 10);
        std::string size_str = std::to_string(quadtree.size());
        DrawText(size_str.c_str(), 10, 30, 10, DARKGRAY);
        DrawText("Left click to draw, right to erase, C to clear", 10, 40, 10, DARKGRAY);

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
                          [&](auto it) { points_to_erase.push_back(*it); });
            for (auto point : points_to_erase) {
                quadtree.erase(point);
            }
            DrawRectangle(min_x, min_y, 2 * EraserRadius, 2 * EraserRadius, BLUE);
        } else if (IsKeyReleased(KEY_C)) {
            quadtree.clear();
        }
        rlImGuiEnd();
        EndDrawing();
    }

    rlImGuiShutdown();
    CloseWindow();
    return 0;
}
