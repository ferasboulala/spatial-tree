# spatial-tree
Header only, no-dependency, dynamic and low latency spatial partitioning data structure.

```C++
#include <spatial-tree.h>

...
  // Basic quadtree
  st::spatial_set<double, 2> quadtree;
  // Basic octree
  st::spatial_set<double, 3> octree;
  // Basic octree bounded by the given bounding box
  st::spatial_set<double, 3> bounded_octree({-10, -10, 10, 10});
  // Quadtree that maps points to std::string
  st::spatial_map<double, std::string, 2> quadtree_map;
  // Spatial tree that maps points to strings in 3D with a maximum leaf size of 256, using 32-bit indices and allowing duplicates
  st::spatial_tree<float, std::string, 3, 256, 32, true> fully_customized_tree;

  // Emplace a point
  quadtree.emplace({0.0, 1.0});
  // Erase a point if it exists
  quadtree.erase({0.0, 1.0});
  // Find a point
  quadtree.find({0.0, 1.0});
  // Find all points in the bounding box
  quadtree.find({-10, -10, 10, 10}, [&](auto it) { /** do something here **/ });
  // Find all points in the bounding sphere
  quadtree.find({5.0, {0.0, 0.0}}, [&](auto it) { /** do something here **/ });
  // Find the 5 nearest points to the provided point
  quadtree.nearest({0.0, 0.0}, 5);
```

## Purpose
- Spatially sort points and optionally store data like in a hash table.
- Need to do nearest neighbour, bounding box range query efficiently.
- Add and remove points without losing query performance or having to reconstruct the entire tree.
- Tune tree parameters (like leaf size and index bit width).

## How to Use
Copy `spatial-tree/spatial-tree.h` in your project and include it.

## Examples
See `example.cpp` and `demo.cpp` for an interactive demo. The demo is also available here: https://ferasboulala.github.io/wasm/draw.html
