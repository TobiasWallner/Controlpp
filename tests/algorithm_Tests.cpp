// google test
#include <gtest/gtest.h>

#include <controlpp.hpp>
#include <Eigen/Dense>

TEST(algorithm, find_enclosing_in_range){
    using namespace controlpp;

    Eigen::VectorXd range(5);
    range << 1.0, 2.0, 3.0, 4.0, 5.0;

    const auto result = find_enclosing(range, 2.3);
    ASSERT_TRUE(result.has_value());
    if(result.has_value()){
        ASSERT_EQ(*result.value().first, 2.0);
        ASSERT_EQ(*result.value().second, 3.0);
    }
}

TEST(algorithm, find_enclosing_before_range){
    using namespace controlpp;

    Eigen::VectorXd range(5);
    range << 1.0, 2.0, 3.0, 4.0, 5.0;

    const auto result = find_enclosing(range, 0.1);
    ASSERT_FALSE(result.has_value());
}

TEST(algorithm, find_enclosing_after_range){
    using namespace controlpp;

    Eigen::VectorXd range(5);
    range << 1.0, 2.0, 3.0, 4.0, 5.0;

    const auto result = find_enclosing(range, 5.1);
    ASSERT_FALSE(result.has_value());
}

TEST(algorithm, shift_up){
    Eigen::Vector<int, 5> v(1, 2, 3, 4, 5);

    controlpp::shift_up(v.data(), v.data() + v.size(), 42);

    ASSERT_EQ(v(0), 42);
    ASSERT_EQ(v(1), 1);
    ASSERT_EQ(v(2), 2);
    ASSERT_EQ(v(3), 3);
    ASSERT_EQ(v(4), 4);
}