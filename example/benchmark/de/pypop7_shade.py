import numpy as np
import time
from pypop7.optimizers.de.shade import SHADE
from pypop7.benchmarks.base_functions import rosenbrock  # Rosenbrock function

# Parameter settings
ndim = 30  # Dimensionality
max_fe = 30000  # Maximum function evaluations
n_runs = 10  # Number of independent runs
population_size = 100  # Population size

# Arrays to store results
best_fitness = np.zeros(n_runs)
run_times = np.zeros(n_runs)

# Optimization problem configuration
problem = {
    'fitness_function': rosenbrock,  # Rosenbrock function
    'ndim_problem': ndim,
    'lower_boundary': -100 * np.ones(ndim),  # Adjusted boundaries
    'upper_boundary': 100 * np.ones(ndim)
}

# Main loop
for i in range(n_runs):
    # Algorithm settings
    options = {
        'max_function_evaluations': max_fe,
        'n_individuals': population_size,
        'seed': np.random.randint(1e5),  # Random seed
        'verbose': False  # Disable internal output
    }
    
    # Run optimization
    start_time = time.time()
    optimizer = SHADE(problem, options)
    results = optimizer.optimize()
    run_times[i] = time.time() - start_time
    
    # Record the best fitness
    best_fitness[i] = results['best_so_far_y']

# Calculate result statistics
mean_fitness = np.mean(best_fitness)
var_fitness = np.var(best_fitness)
mean_time = np.mean(run_times)

# Output results

print(f"Results of SHADE on 30-dimensional Rosenbrock function over {n_runs} runs:")
print(f"Mean best fitness: {mean_fitness:.3e}")
print(f"Variance of best fitness: {var_fitness:.3e}")
print(f"Average runtime: {mean_time:.2f} seconds")