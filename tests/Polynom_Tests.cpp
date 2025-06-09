// google test
#include <gtest/gtest.h>

// controlpp
#include <controlpp/Polynom.hpp>

TEST(Polynom, default_construction){
    control::Polynom<float, 3>();
    control::Polynom<double, 10>();
    control::Polynom<int, 1>();
}

TEST(Polynom, array_construction){
    control::Polynom<int, 3> p({1,2,3});

    ASSERT_EQ(p[0], 1);
    ASSERT_EQ(p[1], 2);
    ASSERT_EQ(p[2], 3);

    ASSERT_EQ(p.order(), 2);
}


TEST(Polynom, vector_construction){
    Eigen::Vector<float, 3> vector({1.0f,2.0f,3.0f});
    control::Polynom<float, 3> p(vector);

    ASSERT_EQ(p.at(0), 1.0f);
    ASSERT_EQ(p.at(1), 2.0f);
    ASSERT_EQ(p.at(2), 3.0f);
}


TEST(Polynom, addition){
    const int arr_a[]  = {1, 2, 3};
    const int arr_b[]  = {4, 5, 6};
    control::Polynom<int, 3> a(arr_a);
    control::Polynom<int, 3> b(arr_b);

    const auto result = a + b;

    for(int i = 0; i <= 2; ++i){
        ASSERT_EQ(result[i], (arr_a[i] + arr_b[i]));
    }
}

TEST(Polynom, scalar_multiplication){
    control::Polynom<int, 3> a({1, 2, 3});
    control::Polynom<int, 3> expected({2, 4, 6});

    const auto result1 = a * 2;
    const auto result2 = 2 * a;
    const auto result3 = a * 3;

    ASSERT_EQ(result1, expected);
    ASSERT_EQ(result2, expected);
    ASSERT_NE(result3, expected);
}

TEST(Polynom, multiplication){
    const control::Polynom<int, 3> a({1, 2, 3});
    const control::Polynom<int, 2> b({1, 2});

    // 1, 2, 3
    //    2, 4, 6
    // --------------
    // 1, 4, 7, 6
    const control::Polynom<int, 4> expected({1, 4, 7, 6});

    const auto result = a * b;

    ASSERT_EQ(result.order(), expected.order());
    ASSERT_EQ(result, expected);
}

