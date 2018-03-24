#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../include/JSO.hpp"
#include "catch.hpp"
#include "integration_tests.hpp"

TEST_CASE("Integration Test with classless fucntion")
{
    srand((unsigned)time(NULL));
    JSO::JSO algorithm(rastrigin_func, 2, -5.12, 5.12);
    algorithm.run();
}

TEST_CASE("Function pointer test with member function")
{
    srand((unsigned)time(NULL));
    TestFunction sphere;
    double f;
    std::vector<double> x = {3.4, 3.0};
    std::function<void(TestFunction*, double*, double*)> fitness_function = &TestFunction::sphere_func;
    fitness_function(&sphere, &x[0], &f);
    REQUIRE(f == 20.56);
}
