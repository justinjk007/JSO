#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../include/jSO.hpp"
#include "catch.hpp"

TEST_CASE("Integration Test")
{
    // random seed is selected based on time according to competition rules
    srand((unsigned)time(NULL));

    jSO::SearchAlgorithm* algorithm = new jSO::LSHADE();

    algorithm->g_problem_size = 2;  // dimension size.
    algorithm->g_max_num_evaluations =
        algorithm->g_problem_size * 10000;  // maximum number of evaluations
    algorithm->g_pop_size =
        (int)round(sqrt(algorithm->g_problem_size) * log(algorithm->g_problem_size) * 25);
    algorithm->g_memory_size = 5;
    algorithm->g_p_best_rate = 0.25;
    algorithm->g_arc_rate    = 1;
    algorithm->domain_min    = -5.12;
    algorithm->domain_max    = 5.12;

    algorithm->run();
    delete algorithm;
}
