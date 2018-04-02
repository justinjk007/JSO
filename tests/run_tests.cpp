#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"
#define private public    // For unit testing private methods
#define protected public  // For unit testing protected methods
#include "../include/JSO.hpp"
#include "integration_tests.hpp"

TEST_CASE("Integration Test with classless fucntion")
{
    srand((unsigned)time(NULL));
    JSO::JSO algorithm(rastrigin_func, 2, -5.12, 5.12);
    // algorithm.run();
}

TEST_CASE("Function pointer test with member function")
{
    srand((unsigned)time(NULL));
    TestFunction sphere;
    double f;
    std::vector<double> x = {3.4, 3.0};
    std::function<void(TestFunction*, double*, double*)> fitness_function =
        &TestFunction::sphere_func;
    fitness_function(&sphere, &x[0], &f);
    REQUIRE(f == 20.56);
}

TEST_CASE("Integration test case with member function")
{
    srand((unsigned)time(NULL));
    TestFunction sphere;
    JSO::JSO algorithm([&sphere](double* p1, double* p2) { sphere.sphere_func(p1, p2); }, 2, -100,
                       100);  // lambda blackmagic
    // algorithm.run();
}

TEST_CASE("Testing randDouble()")
{
    JSO::JSO algorithm(rastrigin_func, 2, -5.12, 5.12);
    double calc[5] = {0};
    for (int i = 0; i < 5; ++i) calc[i] = algorithm.randDouble();
    for (int i = 0; i < 5; ++i) {
        bool cond = calc[i] <= 1 && calc[i] >= 0;
        REQUIRE(cond == true);
    }
}

TEST_CASE("Testing makeNewIndividual()")
{
    JSO::JSO algorithm(rastrigin_func, 2, -5.12, 5.12);
    double* calc[5] = {0};
    for (int i = 0; i < 5; ++i) calc[i] = algorithm.makeNewIndividual();
    for (int i = 0; i < 5; ++i) {
	bool cond = *calc[0] <= 5.12 && *calc[1] >= -5.12;
	REQUIRE(sizeof(calc[i]) == 8);
	REQUIRE(cond == true);
    }
}
