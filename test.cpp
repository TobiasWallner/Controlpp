// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/Polynomial.hpp>


TEST(HalloWorld, hallo_world_test) {
    ASSERT_EQ(true, true);
}


TEST(Polynomial, default_construction){
    control::Polynomial<float, 3>();
    control::Polynomial<double, 10>();
    control::Polynomial<int, 1>();
}

TEST(Polynomial, array_construction){
    control::Polynomial<int, 3> p({1,2,3});

    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);

    ASSERT_EQ(p.order(), 2);
}


TEST(Polynomial, vector_construction){
    Eigen::Vector<float, 3> vector({1.0f,2.0f,3.0f});
    control::Polynomial<float, 3> p(vector);

    ASSERT_EQ(p.at(0), 1.0f);
    ASSERT_EQ(p.at(1), 2.0f);
    ASSERT_EQ(p.at(2), 3.0f);
}


TEST(Polynomial, addition){
    const int arr_a[]  = {1, 2, 3};
    const int arr_b[]  = {4, 5, 6};
    control::Polynomial<int, 3> a(arr_a);
    control::Polynomial<int, 3> b(arr_b);

    const auto result = a + b;

    for(int i = 0; i <= 2; ++i){
        ASSERT_EQ(result[i], (arr_a[i] + arr_b[i]));
    }
}

TEST(Polynomial, scalar_multiplication){
    control::Polynomial<int, 3> a({1, 2, 3});
    control::Polynomial<int, 3> expected({2, 4, 6});

    const auto result1 = a * 2;
    const auto result2 = 2 * a;
    const auto result3 = a * 3;

    ASSERT_EQ(result1, expected);
    ASSERT_EQ(result2, expected);
    ASSERT_NE(result3, expected);
}

TEST(Polynomial, multiplication){
    const control::Polynomial<int, 3> a({1, 2, 3});
    const control::Polynomial<int, 2> b({1, 2});

    // 1, 2, 3
    //    2, 4, 6
    // --------------
    // 1, 4, 7, 6
    const control::Polynomial<int, 4> expected({1, 4, 7, 6});

    const auto result = a * b;

    ASSERT_EQ(result.order(), expected.order());
    ASSERT_EQ(result, expected);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}