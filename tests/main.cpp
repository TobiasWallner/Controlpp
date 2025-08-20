// google test
#include <gtest/gtest.h>

#include <iostream>
#include <exception>
#include <stdexcept>
#include <new>

#include <controlpp.hpp>

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
