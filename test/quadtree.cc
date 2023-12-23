#include "spatial-tree/quadtree.hh"


#include "test/fixture.hh"

#include <gtest/gtest.h>
#include <string>

TEST(QuadTree, Playground) {
    st::QuadTree<const char*, double, 4> tree({-1, 1, 1, -1});
    ASSERT_EQ(tree.size(), 0);
    ASSERT_TRUE(tree.empty());
    tree.clear();
    ASSERT_EQ(tree.size(), 0);
    ASSERT_TRUE(tree.empty());

    const auto ret = tree.emplace(-0.5, 0.5, "point1");
    ASSERT_TRUE(ret.second);
    const auto [x, y, value] = *ret.first;
    ASSERT_STREQ(value, "point1");
    ASSERT_EQ(x, -0.5);
    ASSERT_EQ(y, 0.5);
    ASSERT_FALSE(tree.emplace(-0.5, 0.5, "should not insert").second);
    ASSERT_TRUE(tree.emplace(-0.25, 0.5, "point2").second);
    ASSERT_TRUE(tree.emplace(-0.75, 0.5, "point3").second);
    ASSERT_TRUE(tree.emplace(-0.5, 0.25, "point4").second);
    ASSERT_TRUE(tree.emplace(-0.5, 0.75, "point5").second);

    std::cout << tree.size() << "\n";
}