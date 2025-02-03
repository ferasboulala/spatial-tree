#include <raylib.h>

#include <algorithm>
#include <iostream>
#include <random>

#include "spatial-tree.h"

using CoordinateType = float;

class my_spatial_set : public st::internal::spatial_tree<CoordinateType, Color, 2, 1, 32> {
public:
    void draw(Image& img) const {
        const std::function<void(uint64_t, const st::bounding_box<CoordinateType, 2>&)>
            walk_recursively = [&](auto branch_index, auto boundary) {
                assert(branch_index < this->branches_.size());
                auto& branch = this->branches_[branch_index];

                if (branch.is_terminal()) {
                    return;
                }

                auto [x, y] = boundary.origin();
                ImageDrawLine(&img, x, boundary.starts[1], x, boundary.stops[1], GRAY);
                ImageDrawLine(&img, boundary.starts[0], y, boundary.stops[0], y, GRAY);

                st::internal::unroll_for<BranchingFactor>([&](auto child) {
                    const uint64_t child_index = branch.index_of_first_child + child;
                    const auto     new_boundary = boundary.qrecurse(child);
                    walk_recursively(child_index, new_boundary);
                });
            };

        walk_recursively(0, this->boundary_);
    }
};

static inline Color rainbow(double value) {
    value = std::clamp<double>(value, 0, 1);
    double h = (1 - value) * 240;

    double x = 1 * (1 - fabsf(fmodf(h / 60.0, 2) - 1));
    double r = 0, g = 0, b = 0;
    if (h >= 0 && h < 60) {
        r = 1, g = x, b = 0;
    } else if (h >= 60 && h < 120) {
        r = x, g = 1, b = 0;
    } else if (h >= 120 && h < 180) {
        r = 0, g = 1, b = x;
    } else if (h >= 180 && h < 240) {
        r = 0, g = x, b = 1;
    } else if (h >= 240 && h < 300) {
        r = x, g = 0, b = 1;
    } else if (h >= 300 && h < 360) {
        r = 1, g = 0, b = x;
    }

    Color color;
    color.r = (unsigned char)(r * 255);
    color.g = (unsigned char)(g * 255);
    color.b = (unsigned char)(b * 255);
    color.a = 255;

    return color;
}

int main() {
    static constexpr uint64_t WindowWidth = 800;
    static constexpr uint64_t WindowHeight = 800;
    static constexpr double   BrushRateOfChange = 10;

    my_spatial_set quadtree({0, 0, WindowWidth, WindowHeight});

    std::vector<std::array<CoordinateType, 2> > points_to_erase;
    std::default_random_engine                  gen;
    std::normal_distribution<double>            normal(0.0, 1.0);
    uint64_t                                    brush_size = 20;
    int64_t                                     brush_radius = 20;
    double                                      color_counter = 0;

    SetTraceLogLevel(TraceLogLevel::LOG_FATAL);
    InitWindow(WindowWidth, WindowHeight, "draw");
    SetTargetFPS(60);

    BeginDrawing();
    ClearBackground(RAYWHITE);
    EndDrawing();

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(RAYWHITE);

        Image img = GenImageColor(GetScreenWidth(), GetScreenHeight(), WHITE);
        for (auto [pos, color] : quadtree) {
            ImageDrawPixel(&img, pos[0], pos[1], color);
        }

        brush_size = std::max<uint64_t>(1, M_PI * brush_radius * brush_radius / 1000.0);
        auto mouse_position = GetMousePosition();
        if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
            for (uint64_t i = 0; i < brush_size; ++i) {
                auto offset_x = normal(gen) * brush_radius / 2;
                auto offset_y = normal(gen) * brush_radius / 2;

                std::array<CoordinateType, 2> point = {CoordinateType(mouse_position.x + offset_x),
                                                       CoordinateType(mouse_position.y + offset_y)};
                if (quadtree.fits(point)) {
                    quadtree.emplace(point, rainbow(std::fmod(color_counter, 2) / 2));
                }
            }
            color_counter += GetFrameTime();
        }

        points_to_erase.clear();
        quadtree.find(st::sphere<CoordinateType, 2>{static_cast<CoordinateType>(brush_radius),
                                                    {mouse_position.x, mouse_position.y}},
                      [&](auto it) { points_to_erase.push_back((*it).first); });

        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT)) {
            for (auto point : points_to_erase) {
                quadtree.erase(point);
            }
        } else {
            for (auto point : points_to_erase) {
                ImageDrawPixel(&img, point[0], point[1], RED);
            }
        }

        Vector2 mouse_wheel_movement = GetMouseWheelMoveV();
        if (mouse_wheel_movement.y) {
            brush_radius -= BrushRateOfChange * mouse_wheel_movement.y;
            brush_radius = std::max<int64_t>(1, brush_radius);
        }

        if (IsKeyDown(KEY_D)) {
            quadtree.draw(img);
        }

        Texture tex = LoadTextureFromImage(img);
        UnloadImage(img);
        DrawTexture(tex, 0, 0, WHITE);

        DrawFPS(10, 10);
        const std::string size_str = std::to_string(quadtree.size());
        DrawText(size_str.c_str(), 10, 30, 10, DARKGRAY);
        DrawText("LMB to draw, RMB to erase, wheel to increase brush size", 10, 40, 10, DARKGRAY);

        if (!quadtree.empty()) {
            auto [nearest, _] = *quadtree.nearest({mouse_position.x, mouse_position.y}, 1).front();
            DrawLine(mouse_position.x, mouse_position.y, nearest[0], nearest[1], GRAY);
        }

        DrawCircleLines(mouse_position.x, mouse_position.y, brush_radius, GRAY);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}
