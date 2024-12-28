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
  // Spatial tree that maps points to strings in 3D with a maximum leaf size of 256, using 32-bit indices and checks for duplicates
  st::spatial_tree<float, std::string, 3, 256, 32, true> performance_tuned_tree;

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

  // Range based for-loop
  for (auto& [coordinates, data] : quadtree_map) { /** do something here **/ }
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

## Benchmarks (plugged-in M1 MAX)
Benchmarks with thousands of items per second means once the tree is constructed, that is how many "things" we can do per second. For the rest of the benchmarks like insertions, it refers to how many trees of that size we can construct per second.
To reproduce with different sizes, compile the project using `./runbuild.sh` and then run `./bin/bench`.
| Benchmark                                        | User Counters                               |
|--------------------------------------------------|---------------------------------------------|
| insertions/1048576                               | items_per_second=40.4066/s                  |
| insertions_check_duplicates/1048576              | items_per_second=22.1092/s                  |
| insertions_duplicate/1048576                     | items_per_second=8.67602M/s                 |
| deletions/1048576                                | items_per_second=7.63367/s                  |
| deletions_check_duplicates/1048576               | items_per_second=9.78012/s                  |
| deletions_non_existent/1048576                   | items_per_second=385.227M/s                 |
| deletions_non_existent_check_duplicates/1048576  | items_per_second=300.778M/s                 |
| find/1048576                                     | found=9.467k, items_per_second=69.026k/s    |
| find_sphere/1048576                              | found=30.192k, items_per_second=35.4991k/s  |
| find_single/1048576                              | items_per_second=9.75882M/s                 |
| contains/1048576                                 | items_per_second=178.967M/s                 |
| nearest/1048576                                  | found=16.7772M, items_per_second=2.67842M/s |
| iteration/1048576                                | items_per_second=1.16285G/s                 |
