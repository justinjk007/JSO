#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../include/JSO.hpp"
#include "catch.hpp"
#include "integration_tests.hpp"

TEST_CASE("Integration Test")
{
    // random seed is selected based on time according to competition rules
    srand((unsigned)time(NULL));
    JSO::JSO algorithm(rastrigin_func,2,-5.12,5.12);
    algorithm.run();
}
