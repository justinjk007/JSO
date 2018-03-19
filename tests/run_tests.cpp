#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "../include/jSO.hpp"
#include "catch.hpp"

// Global Variables for DE library
int g_problem_size;
int unsigned g_max_num_evaluations;
int g_pop_size;
double g_arc_rate;
int g_memory_size;
double g_p_best_rate;
double domain_max;
double domain_min;

TEST_CASE("Integration Test")
{
    // random seed is selected based on time according to competition rules
    srand((unsigned)time(NULL));
    g_problem_size        = 2;                       // dimension size.
    g_max_num_evaluations = g_problem_size * 10000;  // available number of fitness evaluations
    g_pop_size            = (int)round(sqrt(g_problem_size) * log(g_problem_size) * 25);
    g_memory_size         = 5;
    g_arc_rate            = 1;
    g_p_best_rate         = 0.25;
    domain_min            = -5.12;
    domain_max            = 5.12;

    jSO::searchAlgorithm* algorithm = new jSO::LSHADE();
    algorithm->run();
    delete algorithm;
}
