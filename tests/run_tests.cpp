#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../include/JSO.hpp"
#include "catch.hpp"

TEST_CASE("Integration Test")
{
    // random seed is selected based on time according to competition rules
    srand((unsigned)time(NULL));

    JSO::JSO algorithm(2,-5.12,5.12);
    algorithm.run();
}
