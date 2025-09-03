// std
#include <numbers>

// google test
#include <gtest/gtest.h>

#include <controlpp.hpp>

TEST(Generator, SineGenerator){
    const float frequency = 10.f;
    const float omega = frequency * 2.f * std::numbers::pi_v<float>;
    const float Tp = 1.f/frequency;
    const float sample_time = Tp/10.f;

    controlpp::SineGenerator sineGen(omega, sample_time);

    for(float time = 0.f; time < Tp; time += sample_time){
        ASSERT_NEAR(sineGen.sin(), std::sin(time * omega), 0.01) << "at time: " << time;
        ASSERT_NEAR(sineGen.cos(), std::cos(time * omega), 0.01) << "at time: " << time;
        sineGen.next();
    }
}