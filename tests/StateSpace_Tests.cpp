// google test
#include <gtest/gtest.h>

// eigen
#include <Eigen/Core>

// controlpp
#include <controlpp/StateSpace.hpp>

TEST(StateSpace, construct_with_to_state_space){
    using namespace controlpp;
    const TransferFunction<float, 1, 2> rp({2.f, 1.f}, {3.f, 2.f, 1.f});
    
    const StateSpace ss = to_state_space(rp);

    const Eigen::Matrix<float, 2, 2> expected_A({
        {0.f, 1.f}, 
        {-3.f, -2.f}
    });
    const Eigen::Matrix<float, 2, 1> expected_B({
        0.f, 
        1.f
    });
    const Eigen::Matrix<float, 1, 2> expected_C({
        2.f, 1.f
    });
    const Eigen::Matrix<float, 1, 1> expected_D(0.f);

    ASSERT_EQ(ss.A(), expected_A);
    ASSERT_EQ(ss.B(), expected_B);
    ASSERT_EQ(ss.C(), expected_C);
    ASSERT_EQ(ss.D(), expected_D);

}