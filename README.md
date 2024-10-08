# spatial-tree
Header only, no-dependency spatial partitioning data structure (a Quadtree in this case) built for low latency dynamic use cases.

## Purpose
- Need to do nearest neighbour, bounding box range query efficiently.
- Add and remove points without losing query performance or having to reconstruct the entire tree (Quadtree invariants are always maintained).

## How to Use
Copy `spatial-tree.h` in your project and include it.

## Examples
See `example.cpp` and `draw.cpp` for an interactive demo.