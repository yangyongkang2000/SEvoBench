#include <iostream>
#include"parameter_observer.h"
#include <fstream>
#include <iostream>
#include <numeric> // For std::accumulate

template <std::floating_point T>
void export_f_cr(const sevobench::experiment::jade_parameter_observer<T>& observer,
                 int problem_index, int instance,
                 const std::string& f_filename, const std::string& cr_filename) {
    // Get the parameter data for the specified problem and instance
    auto instance_data = observer.at(problem_index, instance);

    // Check if there is any data
    if (instance_data.empty()) {
        std::cerr << "No data found for problem_index=" << problem_index
                  << ", instance=" << instance << std::endl;
        return;
    }

    // Determine the number of iterations
    size_t num_iterations = instance_data[0].size(); // Number of iterations per run

    // Prepare vectors to store the sum of f and cr values for each iteration
    std::vector<T> f_sum(num_iterations, T(0));
    std::vector<T> cr_sum(num_iterations, T(0));

    // Sum f and cr values across all runs
    for (const auto& run_data : instance_data) {
        for (size_t i = 0; i < run_data.size(); ++i) {
            f_sum[i] += run_data[i][0]; // Sum f values
            cr_sum[i] += run_data[i][1]; // Sum cr values
        }
    }

    // Calculate the average values
    std::vector<T> f_avg(num_iterations);
    std::vector<T> cr_avg(num_iterations);
    for (size_t i = 0; i < num_iterations; ++i) {
        f_avg[i] = f_sum[i] / instance_data.size(); // Average f value
        cr_avg[i] = cr_sum[i] / instance_data.size(); // Average cr value
    }

    // Open files for writing
    std::ofstream f_file(f_filename);
    std::ofstream cr_file(cr_filename);

    if (!f_file.is_open() || !cr_file.is_open()) {
        std::cerr << "Failed to open output files." << std::endl;
        return;
    }

    // Write f averages to f_filename
    for (const auto& f : f_avg) {
        f_file << f << "\n";
    }

    // Write cr averages to cr_filename
    for (const auto& cr : cr_avg) {
        cr_file << cr << "\n";
    }

    f_file.close();
    cr_file.close();

    std::cout << "Average f values exported to " << f_filename << std::endl;
    std::cout << "Average cr values exported to " << cr_filename << std::endl;
}
int main() {
    auto p= std::make_unique<sevobench::de_module::trackable_jade_parameter<float>>();
    auto m=std::make_unique<sevobench::de_module::ttpb1_mutation<float>>();
    auto c= std::make_unique<sevobench::de_module::binomial_crossover<float>>();
    auto h=std::make_unique<sevobench::de_module::midpoint_target_repair<float>>();
    auto s= std::make_unique<sevobench::de_module::de_population<float>>();
    auto jade_param=p.get();
    auto jade=sevobench::de_module::de_algorithm_builder().
            mutation(std::move(m)).
            parameter(std::move(p)).
            population_strategy(std::move(s)).
            crossover(std::move(c)).
            constraint_handler(std::move(h)).
            build();
    constexpr int dim=30;
    constexpr int pop_size=100;
    constexpr int max_fes=10000*dim;
    constexpr int runs=30;
    auto suite=sevobench::problem::suite_builder<sevobench::problem::cec2017>().dim<30>().instance_count(1).problem_index({11}).build();
    sevobench::experiment::jade_parameter_observer<float> obs(suite,max_fes,30,jade_param,pop_size);
    sevobench::experiment::evo_bench<false>([&](auto &&f){
        sevobench::evolutionary_algorithm alg(max_fes,pop_size,dim);
        sevobench::population<float> pop(pop_size,dim,float(-100),float(100));
        jade_param->reset();
        jade.run(pop,f,f.lower_bound(),f.upper_bound(),alg);
    },suite,obs,runs);
    export_f_cr(obs, 11, 1, "f_values.txt", "cr_values.txt");
    return 0;
}
