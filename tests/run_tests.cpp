#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../src/jSO.hpp"
#include "catch.hpp"
#include "integration_tests.hpp"

TEST_CASE("Sample Test")  // Delete this
{
    int a = 1, b = 0;
    REQUIRE(b == 0);
    REQUIRE(a == Approx(1).epsilon(0.0001));  // Epsilon is tolerance
}
