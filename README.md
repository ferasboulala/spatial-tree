# spatial-tree
Header only, no-dependency spatial partitioning data structure (a Quadtree in this case) built for low latency dynamic use cases.

## Purpose
- Need to do nearest neighbour, bounding box range query efficiently.
- Add and remove points without losing query performance or having to reconstruct the entire tree (Quadtree invariants are always maintained).

## How to Use
Copy `spatial-tree.h` in your project and include it.

## Examples
See `example.cpp` and `draw.cpp` for an interactive demo.

## TODO
- [ ] ND rank (currently 2D only).
- [ ] More optimization low hanging fruits.
- [ ] Add `modify` for moving points which is faster than `emplace` followed by `erase`.
- [ ] Documentation.
- [ ] Share benchmarks.
