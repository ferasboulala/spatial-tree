# spatial-tree
Header only, no-dependency spatial partitioning data structure (a Quadtree in this case) built for low latency dynamic use cases.

## Purpose
- Need to do nearest neighbour, bounding box range query efficiently.
- Add and remove points without having to reconstruct the whole tree.

## How to Use
Copy `spatial-tree.h` in your project and include it.

## Examples
See `example.cpp` and `interactive.cpp`.

## TODO
- [ ] ND rank (currently 2D only).
- [ ] Documentation.