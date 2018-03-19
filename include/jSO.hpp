#ifndef JSO_HPP
#define JSO_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "../tests/integration_tests.hpp"

using namespace std;

#define PI 3.14159265358979323846264338327950288

typedef double variable;
typedef variable* Individual;
typedef double Fitness;

/* extern int g_function_number; */
extern int g_problem_size;
extern unsigned int g_max_num_evaluations;
extern int function_name;

extern int g_pop_size;
extern int g_memory_size;
extern double g_p_best_rate;
extern double g_arc_rate;
extern double domain_max;
extern double domain_min;

extern ofstream outFile;

namespace jSO
{
class searchAlgorithm
{
   public:
    virtual Fitness run() = 0;
    long iteration_number = 0;  // Stores the number of iterations the algorithm has mutated too

   protected:
    void evaluatePopulation(const vector<Individual>& pop, vector<Fitness>& fitness);
    void initializeFitnessFunctionParameters();

    void initializeParameters();
    Individual makeNewIndividual();
    void modifySolutionWithParentMedium(Individual child, Individual parent);
    void setBestSolution(const vector<Individual>& pop, const vector<Fitness>& fitness,
                         Individual& bsf_solution, Fitness& bsf_fitness);

    // Return random value with uniform distribution [0, 1)
    inline double randDouble()
    {
        return (double)rand() / (double)RAND_MAX;
    }

    /*
      Return random value from Cauchy distribution with mean "mu" and variance "gamma"
      http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Cauchy
    */
    inline double cauchy_g(double mu, double gamma)
    {
        return mu + gamma * tan(PI * (randDouble() - 0.5));
    }

    /*
      Return random value from normal distribution with mean "mu" and variance "gamma"
      http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Gauss
    */
    inline double gauss(double mu, double sigma)
    {
        return mu + sigma * sqrt(-2.0 * log(randDouble())) * sin(2.0 * PI * randDouble());
    }

    // Recursive quick sort with index array
    template <class VarType>
    void sortIndexWithQuickSort(VarType array[], int first, int last, int index[])
    {
        VarType x        = array[(first + last) / 2];
        int i            = first;
        int j            = last;
        VarType temp_var = 0;
        int temp_num     = 0;

        while (true) {
            while (array[i] < x) i++;
            while (x < array[j]) j--;
            if (i >= j) break;

            temp_var = array[i];
            array[i] = array[j];
            array[j] = temp_var;

            temp_num = index[i];
            index[i] = index[j];
            index[j] = temp_num;

            i++;
            j--;
        }

        if (first < (i - 1)) sortIndexWithQuickSort(array, first, i - 1, index);
        if ((j + 1) < last) sortIndexWithQuickSort(array, j + 1, last, index);
    }

    int problem_size;
    variable max_region;
    variable min_region;
    Fitness optimum;          // the goal fitness to be reached
    Fitness epsilon;          // acceptable error value
    unsigned int max_num_evaluations;  // max number of evaluations
    int pop_size;             // population size
};

class LSHADE : public searchAlgorithm
{
   public:
    virtual Fitness run();
    void setSHADEParameters();
    void reducePopulationWithSort(vector<Individual>& pop, vector<Fitness>& fitness);
    void operateCurrentToPBest1BinWithArchive(const vector<Individual>& pop, Individual child,
                                              int& target, int& p_best_individual,
                                              variable& scaling_factor, variable& cross_rate,
                                              const vector<Individual>& archive, int& arc_ind_count,
                                              unsigned int nfes);

    int arc_size;
    double arc_rate;
    variable p_best_rate;
    int memory_size;
    int reduction_ind_num;
};

void searchAlgorithm::initializeParameters()
{
    // function_number     = g_function_number;
    problem_size        = g_problem_size;
    max_num_evaluations = g_max_num_evaluations;
    pop_size            = g_pop_size;
    initializeFitnessFunctionParameters();
}

void searchAlgorithm::evaluatePopulation(const vector<Individual>& pop, vector<Fitness>& fitness)
{
    for (int i = 0; i < pop_size; i++) {
        this->iteration_number++;
        rastrigin_func(pop[i], &fitness[i]);  // Call fitness function
    }
}

void searchAlgorithm::initializeFitnessFunctionParameters()
{
    // epsilon is an acceptable error value.
    epsilon    = pow(10.0, -8);
    min_region = domain_min;
    max_region = domain_max;
    optimum    = 0;  // The fitness we need to find, I think
}

// set best solution (bsf_solution) and its fitness value (bsf_fitness) in the initial population
void searchAlgorithm::setBestSolution(const vector<Individual>& pop, const vector<Fitness>& fitness,
                                      Individual& bsf_solution, Fitness& bsf_fitness)
{
    int current_best_individual = 0;

    for (int i = 1; i < pop_size; i++) {
        if (fitness[current_best_individual] > fitness[i]) {
            current_best_individual = i;
        }
    }

    bsf_fitness = fitness[current_best_individual];
    for (int i = 0; i < problem_size; i++) {
        bsf_solution[i] = pop[current_best_individual][i];
    }
}

// make new individual randomly
Individual searchAlgorithm::makeNewIndividual()
{
    Individual individual = (variable*)malloc(sizeof(variable) * problem_size);

    for (int i = 0; i < problem_size; i++) {
        individual[i] = ((max_region - min_region) * randDouble()) + min_region;
    }

    return individual;
}

/*
  For each dimension j, if the mutant vector element v_j is outside the boundaries [x_min , x_max],
  we applied this bound handling method
  If you'd like to know that precisely, please read:
  J. Zhang and A. C. Sanderson, "JADE: Adaptive differential evolution with optional external
  archive,"
  IEEE Tran. Evol. Comput., vol. 13, no. 5, pp. 945–958, 2009.
 */
void searchAlgorithm::modifySolutionWithParentMedium(Individual child, Individual parent)
{
    int l_problem_size    = problem_size;
    variable l_min_region = min_region;
    variable l_max_region = max_region;

    for (int j = 0; j < l_problem_size; j++) {
        if (child[j] < l_min_region) {
            child[j] = (l_min_region + parent[j]) / 2.0;
        } else if (child[j] > l_max_region) {
            child[j] = (l_max_region + parent[j]) / 2.0;
        }
    }
}
Fitness LSHADE::run()
{
    cout << scientific << setprecision(8);
    initializeParameters();
    setSHADEParameters();

    // cout << pop_size << endl;
    // cout << arc_size << endl;
    // cout << p_best_rate << endl;
    // cout << memory_size << endl;

    vector<Individual> pop;
    vector<Fitness> fitness(pop_size, 0);
    vector<Individual> children;
    vector<Fitness> children_fitness(pop_size, 0);

    // initialize population
    for (int i = 0; i < pop_size; i++) {
        pop.push_back(makeNewIndividual());
        children.push_back((variable*)malloc(sizeof(variable) * problem_size));
    }

    // evaluate the initial population's fitness values
    evaluatePopulation(pop, fitness);

    Individual bsf_solution = (variable*)malloc(sizeof(variable) * problem_size);
    Fitness bsf_fitness;
    unsigned int nfes = 0;

    if ((fitness[0] - optimum) < epsilon) fitness[0] = optimum;
    bsf_fitness = fitness[0];
    for (int j = 0; j < problem_size; j++) bsf_solution[j] = pop[0][j];
    /////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < pop_size; i++) {
        nfes++;

        if ((fitness[i] - optimum) < epsilon) fitness[i] = optimum;

        if (fitness[i] < bsf_fitness) {
            bsf_fitness = fitness[i];
            for (int j = 0; j < problem_size; j++) bsf_solution[j] = pop[i][j];
        }

        // if (nfes % 1000 == 0) {
        //   //      cout << nfes << " " << bsf_fitness - optimum << endl;
        //   cout << bsf_fitness - optimum << endl;
        // }

        if (nfes >= max_num_evaluations) break;
    }
    ////////////////////////////////////////////////////////////////////////////

    // for external archive
    int arc_ind_count = 0;
    int random_selected_arc_ind;
    vector<Individual> archive;
    for (int i = 0; i < arc_size; i++)
        archive.push_back((variable*)malloc(sizeof(variable) * problem_size));

    int num_success_params     = 0;
    int old_num_success_params = 0;
    vector<variable> success_sf;
    vector<variable> success_cr;
    vector<variable> dif_fitness;

    // the contents of M_f and M_cr are all initialiezed 0.5
    vector<variable> memory_sf(memory_size, 0.3);  // jSO
    // vector <variable> memory_cr(memory_size, 0.5);          // iL-SHADE
    vector<variable> memory_cr(memory_size, 0.8);

    variable temp_sum_sf;
    variable temp_sum_cr;
    variable sum;
    variable weight;

    // memory index counter
    int memory_pos = 0;

    // for new parameters sampling
    variable mu_sf, mu_cr;
    int random_selected_period;
    variable* pop_sf = (variable*)malloc(sizeof(variable) * pop_size);
    variable* pop_cr = (variable*)malloc(sizeof(variable) * pop_size);

    // for current-to-pbest/1
    int p_best_ind;
    int p_num         = round(pop_size * p_best_rate);
    int* sorted_array = (int*)malloc(sizeof(int) * pop_size);
    Fitness* temp_fit = (Fitness*)malloc(sizeof(Fitness) * pop_size);

    // for linear population size reduction
    int max_pop_size = pop_size;
    int min_pop_size = 4;
    int plan_pop_size;

    // main loop
    while (nfes < max_num_evaluations) {
        for (int i = 0; i < pop_size; i++) sorted_array[i] = i;
        for (int i = 0; i < pop_size; i++) temp_fit[i] = fitness[i];
        sortIndexWithQuickSort(&temp_fit[0], 0, pop_size - 1, sorted_array);

        for (int target = 0; target < pop_size; target++) {
            // In each generation, CR_i and F_i used by each individual x_i are generated by first
            // selecting an index r_i randomly from [1, H]
            random_selected_period = rand() % memory_size;
            if (random_selected_period == memory_size - 1) {
                // mu_sf = 0.9; mu_cr = 0.9;
                mu_sf = 0.9;
                mu_cr = 0.9;
            } else {
                mu_sf = memory_sf[random_selected_period];
                mu_cr = memory_cr[random_selected_period];
            }

            // generate CR_i and repair its value
            // if (mu_cr == -1) {
            if (mu_cr < 0) {           // JANEZ
                pop_cr[target] = 0.0;  // LSHADE 0
            } else {
                pop_cr[target] = gauss(mu_cr, 0.1);
                if (pop_cr[target] > 1)
                    pop_cr[target] = 1;
                else if (pop_cr[target] < 0)
                    pop_cr[target] = 0;
            }
            //   if (nfes< 0.25*max_num_evaluations && pop_cr[target] < 0.5) pop_cr[target] = 0.5;
            //   // iL-SHADE
            //   if (nfes< 0.50*max_num_evaluations && pop_cr[target] < 0.25) pop_cr[target] = 0.25;
            //   // iL-SHADE
            if (nfes < 0.25 * max_num_evaluations && pop_cr[target] < 0.7)
                pop_cr[target] = 0.7;  // jSO
            if (nfes < 0.50 * max_num_evaluations && pop_cr[target] < 0.6)
                pop_cr[target] = 0.6;  // jSO

            // generate F_i and repair its value
            do {
                pop_sf[target] = cauchy_g(mu_sf, 0.1);
            } while (pop_sf[target] <= 0.0);

            if (pop_sf[target] > 1) pop_sf[target] = 1.0;
            //  if (nfes< 0.25*max_num_evaluations && pop_sf[target] > 0.8) pop_sf[target] = 0.8;
            //  // iL-SHADE
            //  if (nfes< 0.50*max_num_evaluations && pop_sf[target] > 0.9) pop_sf[target] = 0.9;
            //  // iL-SHADE
            if (nfes < 0.6 * max_num_evaluations && pop_sf[target] > 0.7)
                pop_sf[target] = 0.7;  // jSO

            // //p-best individual is randomly selected from the top pop_size *  p_i members
            //      p_best_ind = sorted_array[rand() % p_num]; // L-SHADE
            // p-best individual is randomly selected from the top pop_size *  p_i members
            do {
                p_best_ind = sorted_array[rand() % p_num];
            } while (nfes < 0.50 * max_num_evaluations && p_best_ind == target);  // iL-SHADE

            operateCurrentToPBest1BinWithArchive(pop, &children[target][0], target, p_best_ind,
                                                 pop_sf[target], pop_cr[target], archive,
                                                 arc_ind_count, nfes);

            // print info
            // cout << nfes <<" "<< target <<" "<< pop_sf[target] <<" "<< pop_cr[target] << endl;
        }

        // evaluate the children's fitness values
        evaluatePopulation(children, children_fitness);

        /////////////////////////////////////////////////////////////////////////
        // update the bsf-solution and check the current number of fitness evaluations
        // if the current number of fitness evaluations over the max number of fitness evaluations,
        // the search is terminated
        // So, this program is unconcerned about L-SHADE algorithm directly
        for (int i = 0; i < pop_size; i++) {
            nfes++;

            // following the rules of CEC 2014 real parameter competition,
            // if the gap between the error values of the best solution found and the optimal
            // solution was 10^{−8} or smaller,
            // the error was treated as 0
            if ((children_fitness[i] - optimum) < epsilon) children_fitness[i] = optimum;

            if (children_fitness[i] < bsf_fitness) {
                bsf_fitness = children_fitness[i];
                for (int j = 0; j < problem_size; j++) bsf_solution[j] = children[i][j];
            }

            // if (nfes % 1000 == 0) {
            // //      cout << nfes << " " << bsf_fitness - optimum << endl;
            // 	cout << bsf_fitness - optimum << endl;
            // }

            // CEC 2014: Record function error value (Fi(x)-Fi(x*)) after (0.01, 0.02, 0.03, 0.05,
            // 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)*MaxFES for each run.
            if (nfes % (g_problem_size * 100) == 0) {
                if (nfes == (int)(0.01 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.02 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.03 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.05 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.1 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.2 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.3 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.4 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.5 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.6 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.7 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.8 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(0.9 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << " ";
                else if (nfes == (int)(1.0 * max_num_evaluations))
                    outFile << bsf_fitness - optimum << endl;
            }

            if (nfes >= max_num_evaluations) break;
        }
        ////////////////////////////////////////////////////////////////////////////

        // generation alternation
        for (int i = 0; i < pop_size; i++) {
            if (children_fitness[i] == fitness[i]) {
                fitness[i] = children_fitness[i];
                for (int j = 0; j < problem_size; j++) pop[i][j] = children[i][j];
            } else if (children_fitness[i] < fitness[i]) {
                dif_fitness.push_back(fabs(fitness[i] - children_fitness[i]));
                fitness[i] = children_fitness[i];
                //  	for (int j = 0; j < problem_size; j ++) pop[i][j] = children[i][j];   //
                //  iL-SHADE

                // successful parameters are preserved in S_F and S_CR
                success_sf.push_back(pop_sf[i]);
                success_cr.push_back(pop_cr[i]);

                // parent vectors x_i which were worse than the trial vectors u_i are preserved
                if (arc_size > 1) {
                    if (arc_ind_count < arc_size) {
                        for (int j = 0; j < problem_size; j++)
                            archive[arc_ind_count][j] = pop[i][j];
                        arc_ind_count++;
                    }
                    // Whenever the size of the archive exceeds, randomly selected elements are
                    // deleted to make space for the newly inserted elements
                    else {
                        random_selected_arc_ind = rand() % arc_size;
                        for (int j = 0; j < problem_size; j++)
                            archive[random_selected_arc_ind][j] = pop[i][j];
                    }
                }
                for (int j = 0; j < problem_size; j++) pop[i][j] = children[i][j];  // jSO
            }
        }

        old_num_success_params = num_success_params;
        num_success_params     = success_sf.size();

        // if numeber of successful parameters > 0, historical memories are updated
        if (num_success_params > 0) {
            variable old_sf, old_cr;  // Janez
            old_sf = memory_sf[memory_pos];
            old_cr = memory_cr[memory_pos];

            memory_sf[memory_pos] = 0;
            memory_cr[memory_pos] = 0;
            temp_sum_sf           = 0;
            temp_sum_cr           = 0;
            sum                   = 0;

            for (int i = 0; i < num_success_params; i++) sum += dif_fitness[i];

            // weighted lehmer mean
            for (int i = 0; i < num_success_params; i++) {
                weight = dif_fitness[i] / sum;

                memory_sf[memory_pos] += weight * success_sf[i] * success_sf[i];
                temp_sum_sf += weight * success_sf[i];

                memory_cr[memory_pos] += weight * success_cr[i] * success_cr[i];
                temp_sum_cr += weight * success_cr[i];
            }

            memory_sf[memory_pos] /= temp_sum_sf;

            if (temp_sum_cr == 0 || memory_cr[memory_pos] == -1)
                memory_cr[memory_pos] = -1;
            else
                memory_cr[memory_pos] /= temp_sum_cr;

            // JANEZ
            memory_sf[memory_pos] = (memory_sf[memory_pos] + old_sf) / 2.0;
            memory_cr[memory_pos] = (memory_cr[memory_pos] + old_cr) / 2.0;

            // if (memory_cr[memory_pos] < 0) {
            //   cout << " nfes: " << nfes ;
            //   cout << " memory_sf[memory_pos] = "<< memory_sf[memory_pos] << "
            //   memory_cr[memory_pos] = " << memory_cr[memory_pos];
            //   cout << " old_sf = "<< old_sf << " old_cr = " << old_cr<< endl;
            //}
            // increment the counter
            memory_pos++;
            if (memory_pos >= memory_size) memory_pos = 0;

            // clear out the S_F, S_CR and delta fitness
            success_sf.clear();
            success_cr.clear();
            dif_fitness.clear();
        }

        // calculate the population size in the next generation
        plan_pop_size = round(
            (((min_pop_size - max_pop_size) / (double)max_num_evaluations) * nfes) + max_pop_size);

        if (pop_size > plan_pop_size) {
            reduction_ind_num = pop_size - plan_pop_size;
            if (pop_size - reduction_ind_num < min_pop_size)
                reduction_ind_num = pop_size - min_pop_size;

            reducePopulationWithSort(pop, fitness);

            // resize the archive size
            arc_size = pop_size * g_arc_rate;
            if (arc_ind_count > arc_size) arc_ind_count = arc_size;

            // resize the number of p-best individuals
            p_best_rate =
                g_p_best_rate * (1.0 - 0.5 * nfes / (double)max_num_evaluations);  // JANEZ
            p_num = round(pop_size * p_best_rate);
            if (p_num <= 1) p_num = 2;
        }
    }

    free(pop_sf);
    free(pop_cr);
    for (int i = 0; i < pop_size; i++) {
        free(pop[i]);       //.clear();  // JANEZ
        free(children[i]);  //.clear();  // JANEZ
    }
    pop.clear();       // JANEZ
    children.clear();  // JANEZ
    fitness.clear();
    children_fitness.clear();
    free(bsf_solution);
    for (int i = 0; i < arc_size; i++) free(archive[i]);
    archive.clear();
    memory_sf.clear();
    // sorted_array.clear();
    free(sorted_array);
    free(temp_fit);  //.clear();

    return bsf_fitness - optimum;
}

void LSHADE::operateCurrentToPBest1BinWithArchive(const vector<Individual>& pop, Individual child,
                                                  int& target, int& p_best_individual,
                                                  variable& scaling_factor, variable& cross_rate,
                                                  const vector<Individual>& archive,
                                                  int& arc_ind_count, unsigned int nfes)
{
    int r1, r2;

    variable jF = scaling_factor;  // jSO
    if (nfes < 0.2 * max_num_evaluations)
        jF = jF * 0.7;  // jSO
    else if (nfes < 0.4 * max_num_evaluations)
        jF = jF * 0.8;  // jSO
    else
        jF = jF * 1.2;  // jSO

    do {
        r1 = rand() % pop_size;
    } while (r1 == target);
    do {
        r2 = rand() % (pop_size + arc_ind_count);
    } while ((r2 == target) || (r2 == r1));

    int random_variable = rand() % problem_size;

    if (r2 >= pop_size) {
        r2 -= pop_size;
        for (int i = 0; i < problem_size; i++) {
            if ((randDouble() < cross_rate) || (i == random_variable)) {
                child[i] = pop[target][i] + jF * (pop[p_best_individual][i] - pop[target][i]) +
                           scaling_factor * (pop[r1][i] - archive[r2][i]);  // jSO
            } else {
                child[i] = pop[target][i];
            }
        }
    } else {
        for (int i = 0; i < problem_size; i++) {
            if ((randDouble() < cross_rate) || (i == random_variable)) {
                child[i] = pop[target][i] + jF * (pop[p_best_individual][i] - pop[target][i]) +
                           scaling_factor * (pop[r1][i] - pop[r2][i]);  // jSO
            } else {
                child[i] = pop[target][i];
            }
        }
    }

    // If the mutant vector violates bounds, the bound handling method is applied
    modifySolutionWithParentMedium(child, pop[target]);
}

void LSHADE::reducePopulationWithSort(vector<Individual>& pop, vector<Fitness>& fitness)
{
    int worst_ind;

    for (int i = 0; i < reduction_ind_num; i++) {
        worst_ind = 0;
        for (int j = 1; j < pop_size; j++) {
            if (fitness[j] > fitness[worst_ind]) worst_ind = j;
        }

        pop.erase(pop.begin() + worst_ind);
        fitness.erase(fitness.begin() + worst_ind);
        pop_size--;
    }
}

void LSHADE::setSHADEParameters()
{
    arc_rate    = g_arc_rate;
    arc_size    = (int)round(pop_size * arc_rate);
    p_best_rate = g_p_best_rate;
    memory_size = g_memory_size;
}

}  // End of namspace JSO
#endif /* JSO_HPP */
